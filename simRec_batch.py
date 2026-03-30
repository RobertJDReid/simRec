"""
simRec_batch.py

Batch runner for simRec. Simulates many independent cell lineages and writes
classified recombination events in tab-separated format.

Each cell is one independent call to run_simulation() followed by
classify_events(). Cells share no state. Output is identical in format to
simRec's --events flag, with an additional 'cell' column.

Output columns:
    cell                : integer cell identifier (1-based)
    time                : ISO timestamp of the batch run start (same for all rows)
    event               : GC-NCO | GC-CO | CO-terminal
    start               : first bp of LOH block
    end                 : last bp of LOH block
    haplotype           : A or B
    left                : LOH state of left neighbor (A, B, het, tel)
    right               : LOH state of right neighbor (A, B, het, tel)
    adjacent_to_terminal: True or False

Usage
-----
    # Write to stdout (default):
    python simRec_batch.py genome_chrII.csv --n-cells 1000 --n-gen 5000

    # Write to file:
    python simRec_batch.py genome_chrII.csv --n-cells 1000 --n-gen 5000 --observed-out results.tsv

    # Pipe into another script:
    python simRec_batch.py genome_chrII.csv --n-cells 500 --n-gen 5000 | downstream.py

Architecture note
-----------------
The per-cell function run_one_cell() is a plain module-level function so it can
be passed to multiprocessing.Pool.map() without pickling issues when
parallelisation is added in the future.
"""

import sys
import random
import argparse
from datetime import datetime

from tqdm import tqdm

from simRec import (
    load_genome,
    run_simulation,
    classify_events,
    reclassify_terminal_clusters,
)

__version__ = "0.10"

# ---------------------------------------------------------------------------
# Per-cell worker
# ---------------------------------------------------------------------------

def run_one_cell(cell_id, genome, n_gen, gc_min, gc_max, co_prob, seed):
    """
    Simulate one cell lineage and return a list of event row tuples.

    Parameters
    ----------
    cell_id : int
        1-based integer identifier for this cell.
    genome : dict
        Pre-replication cell dict {"A": chrom, "B": chrom} from load_genome().
        Each call receives its own deep-copied genome so cells are independent.
    n_gen : int
        Number of generations to simulate.
    gc_min, gc_max : int
        Gene conversion tract size bounds in bp.
    co_prob : float
        Crossover probability per recombination event.
    seed : int
        Per-cell random seed derived from the master seed.

    Returns
    -------
    list of tuple
        One tuple per classified event, ready to join as a TSV row.
    """
    import copy
    rng = random.Random(seed)

    # Give this cell its own RNG state without disturbing the global state
    # by temporarily swapping the module-level random state
    saved_state = random.getstate()
    random.setstate(rng.getstate())

    try:
        cell_genome = copy.deepcopy(genome)
        final_cell, event_log, _ = run_simulation(
            cell_genome,
            n_gen    = n_gen,
            gc_min   = gc_min,
            gc_max   = gc_max,
            co_prob  = co_prob,
        )
        events = reclassify_terminal_clusters(
            classify_events(final_cell, gc_max=gc_max),
            gc_max=gc_max,
        )
    finally:
        random.setstate(saved_state)

    rows = []
    for e in events:
        rows.append((
            str(cell_id),
            e["type"],
            str(e["start"]),
            str(e["end"]),
            e["haplotype"],
            e["flanking_left"],
            e["flanking_right"],
            str(e.get("adjacent_to_terminal", False)),
            str(e.get("reclassified", False)),
            str(e.get("complex", False)),
        ))

    obs_rows = []
    for e in event_log:
        obs_rows.append((
            str(cell_id),
            str(e["gen"]),
            e["chrom"],
            e["type"],
            str(e["start"]),
            str(e["end"]),
        ))

    return rows, obs_rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Batch runner for simRec — simulates many independent cell lineages.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("genome_file",
                        help="CSV with columns: chromosome, length, CEN")
    parser.add_argument("--version", action="version",
                        version=f"simRec_batch {__version__}")
    parser.add_argument("--n-cells", type=int, default=1000, dest="n_cells",
                        help="Number of independent cell lineages to simulate")
    parser.add_argument("--n-gen",   type=int, default=5000, dest="n_gen",
                        help="Number of generations per cell")
    parser.add_argument("--p-rec",   type=float, default=3.7e-10, dest="p_rec",
                        help="Recombination probability per bp per generation "
                             "(default: 3.7e-10, lower-bound estimate from "
                             "Sui et al. 2020 genome-wide MA study)")
    parser.add_argument("--gc-min",  type=int, default=100, dest="gc_min",
                        help="Minimum gene conversion tract size (bp)")
    parser.add_argument("--gc-max",  type=int, default=5000, dest="gc_max",
                        help="Maximum gene conversion tract size (bp)")
    parser.add_argument("--co-nco",  type=float, default=0.5, dest="co_prob",
                        help="Crossover probability per recombination event")
    parser.add_argument("--seed",    type=int, default=None,
                        help="Master random seed for reproducibility")
    parser.add_argument("--observed-out", type=str, default=None, dest="observed_out",
                        help="Output file for post-hoc classified events (observed from the final "
                             "LOH map). If not given, writes to stdout.")
    parser.add_argument("--logged-out", type=str, default=None, dest="logged_out",
                        help="Output file for logged recombination events tracked as the "
                             "simulation runs. Written as TSV with a 'cell' column. "
                             "Columns: cell, gen, chrom, type, start, end")
    args = parser.parse_args()

    # -----------------------------------------------------------------------
    # Set up output streams
    # -----------------------------------------------------------------------
    if args.observed_out:
        observed_fh = open(args.observed_out, "w")
    else:
        observed_fh = sys.stdout

    if args.logged_out:
        logged_fh = open(args.logged_out, "w")
        logged_fh.write("\t".join(("cell", "gen", "chrom", "type", "start", "end")) + "\n")
    else:
        logged_fh = None

    # -----------------------------------------------------------------------
    # Master seed → per-cell seeds
    # -----------------------------------------------------------------------
    master_rng = random.Random(args.seed)
    cell_seeds = [master_rng.randint(0, 2**31 - 1) for _ in range(args.n_cells)]

    # -----------------------------------------------------------------------
    # Load genome once; each cell gets a deep copy inside run_one_cell()
    # -----------------------------------------------------------------------
    genome = load_genome(args.genome_file, p_rec_default=args.p_rec)
    chrom_name   = genome["A"]["name"]
    chrom_length = genome["A"]["length"]

    # -----------------------------------------------------------------------
    # Progress bar goes to stderr so it doesn't pollute stdout/file output
    # -----------------------------------------------------------------------
    timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")

    header = ("cell", "time", "event", "start", "end",
              "haplotype", "left", "right", "adjacent_to_terminal",
              "reclassified", "complex")
    observed_fh.write("\t".join(header) + "\n")

    desc = f"simRec_batch  chr={chrom_name}  L={chrom_length:,}  " \
           f"n_gen={args.n_gen:,}  p_rec={args.p_rec:.2e}"

    try:
        with tqdm(total=args.n_cells, desc=desc, file=sys.stderr,
                  unit="cell", dynamic_ncols=True,
                  disable=args.observed_out is None) as pbar:
            for cell_id in range(1, args.n_cells + 1):
                rows, log_rows = run_one_cell(
                    cell_id  = cell_id,
                    genome   = genome,
                    n_gen    = args.n_gen,
                    gc_min   = args.gc_min,
                    gc_max   = args.gc_max,
                    co_prob  = args.co_prob,
                    seed     = cell_seeds[cell_id - 1],
                )
                for row in rows:
                    observed_fh.write(row[0] + "\t" + timestamp + "\t" + "\t".join(row[1:]) + "\n")
                if logged_fh is not None:
                    for log_row in log_rows:
                        logged_fh.write("\t".join(log_row) + "\n")
                pbar.update(1)
    finally:
        if args.observed_out:
            observed_fh.close()
        if logged_fh is not None:
            logged_fh.close()


if __name__ == "__main__":
    main()
