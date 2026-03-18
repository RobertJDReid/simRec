"""
recombination_sim/simRec.py

Mitotic recombination simulation for a diploid yeast chromosome.

Chromosome slots are represented as lists of (start, end, value) tuples.
Helper functions handle all interval operations.

Cell structure (post-replication):
    cell = {
        "A": [chrA1, chrA2],   # sister chromatid pair from homolog A
        "B": [chrB1, chrB2]    # sister chromatid pair from homolog B
    }

Pre-replication / post-segregation cell:
    cell = {
        "A": chr_A,
        "B": chr_B
    }
"""

import random
import math
import csv
import copy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ---------------------------------------------------------------------------
# Chromosome object
# ---------------------------------------------------------------------------

def make_chromosome(name, length, cen_start, cen_end, p_rec_default=1e-5, parent="A"):
    """
    Create a chromosome object as a dict of named tuple-list slots.

    Slots:
        position  : [(1, length, None)]          — fixed, scalar metadata
        cen       : [(cen_start, cen_end, None)] — fixed protected interval
        p_rec     : [(1, length, p_rec_default)] — recombination probability density
        haplotype : [(1, length, parent)]        — parental identity label
    """
    return {
        "name": name,
        "length": length,
        "position": [(1, length, None)],
        "cen": [(cen_start, cen_end, None)],
        "p_rec": [(1, length, p_rec_default)],
        "haplotype": [(1, length, parent)],
    }


def clone_chromosome(chrom):
    """Deep copy a chromosome object."""
    return copy.deepcopy(chrom)


# ---------------------------------------------------------------------------
# Tuple-list interval helpers
# ---------------------------------------------------------------------------

def query_value(segments, pos):
    """Return the value of the segment containing pos (1-based, inclusive)."""
    for start, end, val in segments:
        if start <= pos <= end:
            return val
    raise ValueError(f"Position {pos} not covered by any segment.")


def replace_interval(segments, new_start, new_end, new_val):
    """
    Replace the value in [new_start, new_end] with new_val.
    Segments outside this interval are preserved.
    Returns a new sorted, merged segment list.
    """
    result = []
    for start, end, val in segments:
        if end < new_start or start > new_end:
            # entirely outside replacement zone
            result.append((start, end, val))
        else:
            # trim left remnant
            if start < new_start:
                result.append((start, new_start - 1, val))
            # trim right remnant
            if end > new_end:
                result.append((new_end + 1, end, val))
    result.append((new_start, new_end, new_val))
    result.sort(key=lambda x: x[0])
    return merge_segments(result)


def merge_segments(segments):
    """Merge adjacent segments that share the same value."""
    if not segments:
        return []
    merged = [segments[0]]
    for start, end, val in segments[1:]:
        prev_start, prev_end, prev_val = merged[-1]
        if val == prev_val and start == prev_end + 1:
            merged[-1] = (prev_start, end, prev_val)
        else:
            merged.append((start, end, val))
    return merged


def total_length(segments):
    return sum(end - start + 1 for start, end, _ in segments)


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def load_genome(filepath, p_rec_default=1e-5, cen_width=200):
    """
    Parse a CSV with columns: chromosome, length, CEN
    CEN column gives the centromere midpoint; cen_width sets the protected interval.

    Returns a pre-replication cell dict:
        { "A": chr_A, "B": chr_B }
    for the first chromosome in the file (single-chromosome mode).
    """
    with open(filepath) as fh:
        reader = csv.DictReader(fh, skipinitialspace=True)
        rows = list(reader)

    if not rows:
        raise ValueError("Genome file is empty.")

    # Single-chromosome mode: use first row only
    row = rows[0]
    name = row["chromosome"].strip()
    length = int(row["length"])
    cen_mid = int(row["CEN"])
    cen_start = max(1, cen_mid - cen_width // 2)
    cen_end = min(length, cen_mid + cen_width // 2)

    chr_A = make_chromosome(name, length, cen_start, cen_end, p_rec_default, parent="A")
    chr_B = make_chromosome(name, length, cen_start, cen_end, p_rec_default, parent="B")

    return {"A": chr_A, "B": chr_B}


# ---------------------------------------------------------------------------
# Replication
# ---------------------------------------------------------------------------

def replicate_cell(pre_rep_cell):
    """
    Duplicate each homolog to produce two sister chromatid pairs.

    Input:  {"A": chr_A, "B": chr_B}
    Output: {"A": [chrA1, chrA2], "B": [chrB1, chrB2]}
    """
    return {
        "A": [clone_chromosome(pre_rep_cell["A"]),
              clone_chromosome(pre_rep_cell["A"])],
        "B": [clone_chromosome(pre_rep_cell["B"]),
              clone_chromosome(pre_rep_cell["B"])],
    }


# ---------------------------------------------------------------------------
# Recombination — helpers
# ---------------------------------------------------------------------------

def _cen_interval(chrom):
    """Return (cen_start, cen_end) from the chromosome's CEN slot."""
    seg = chrom["cen"]
    return seg[0][0], seg[0][1]


def _sample_gc_site(chrom, gc_min=100, gc_max=5000):
    """
    Sample a valid gene conversion tract on chrom.
    - Centered on a uniformly random non-CEN position weighted by p_rec.
    - Clipped to chromosome ends.
    - Rejected if the tract overlaps the centromere.

    Returns (gc_start, gc_end) or None if no valid site found after max_tries.
    """
    length = chrom["length"]
    cen_start, cen_end = _cen_interval(chrom)
    p_rec_segs = chrom["p_rec"]
    max_tries = 1000

    # Build list of sub-intervals that are entirely outside the CEN
    eligible = []
    for start, end, val in p_rec_segs:
        if end < cen_start or start > cen_end:
            # Entirely outside CEN
            eligible.append((start, end, val))
        else:
            # Partially overlaps CEN — keep the non-CEN flanks
            if start < cen_start:
                eligible.append((start, cen_start - 1, val))
            if end > cen_end:
                eligible.append((cen_end + 1, end, val))
    if not eligible:
        return None

    weights = [(end - start + 1) * val for start, end, val in eligible]
    total_weight = sum(weights)

    for _ in range(max_tries):
        # Pick a segment proportional to weight
        r = random.uniform(0, total_weight)
        cumulative = 0
        chosen_seg = eligible[-1]
        for seg, w in zip(eligible, weights):
            cumulative += w
            if r <= cumulative:
                chosen_seg = seg
                break

        center = random.randint(chosen_seg[0], chosen_seg[1])
        gc_size = random.randint(gc_min, gc_max)
        half = gc_size // 2
        gc_start = max(1, center - half)
        gc_end = min(length, gc_start + gc_size - 1)

        # Reject only if the final clipped tract overlaps centromere
        if gc_end >= cen_start and gc_start <= cen_end:
            # Nudge the tract away from the centromere if possible
            if center < cen_start:
                gc_end = min(cen_start - 1, gc_end)
            else:
                gc_start = max(cen_end + 1, gc_start)
            # After nudge, verify tract still has positive length
            if gc_start > gc_end:
                continue

        return gc_start, gc_end

    return None  # failed to place a valid GC (should be extremely rare)


def _apply_gene_conversion(initiator, donor, gc_start, gc_end):
    """
    Copy haplotype from donor [gc_start, gc_end] into initiator.
    Modifies initiator in place.
    """
    donor_val = query_value(donor["haplotype"], (gc_start + gc_end) // 2)
    initiator["haplotype"] = replace_interval(
        initiator["haplotype"], gc_start, gc_end, donor_val
    )


def _apply_crossover(chrom_a, chrom_b, gc_start, gc_end):
    """
    Swap the telomere-distal arm of chrom_a and chrom_b relative to the GC tract.
    - If GC is on the left arm (gc_end < cen_start): swap 1 .. gc_start-1
    - If GC is on the right arm (gc_start > cen_end): swap gc_end+1 .. length
    Modifies both chromosomes in place.
    """
    cen_start, cen_end = _cen_interval(chrom_a)
    length = chrom_a["length"]

    if gc_end < cen_start:
        # Left arm: swap from 1 to gc_start - 1
        swap_start, swap_end = 1, gc_start - 1
    else:
        # Right arm: swap from gc_end + 1 to end
        swap_start, swap_end = gc_end + 1, length

    if swap_start > swap_end:
        return  # degenerate case: GC reaches the telomere, nothing to swap

    # Snapshot the haplotype segments in the swap region for both chromosomes
    def extract_region(chrom, s, e):
        result = []
        for seg_s, seg_e, val in chrom["haplotype"]:
            lo = max(seg_s, s)
            hi = min(seg_e, e)
            if lo <= hi:
                result.append((lo, hi, val))
        return result

    region_a = extract_region(chrom_a, swap_start, swap_end)
    region_b = extract_region(chrom_b, swap_start, swap_end)

    def apply_region(chrom, region):
        for seg_s, seg_e, val in region:
            chrom["haplotype"] = replace_interval(chrom["haplotype"], seg_s, seg_e, val)

    apply_region(chrom_a, region_b)
    apply_region(chrom_b, region_a)


# ---------------------------------------------------------------------------
# Recombination — main function
# ---------------------------------------------------------------------------

def recombine(post_rep_cell, gc_min=100, gc_max=5000):
    """
    Simulate recombination on a post-replication cell.

    For each chromosome:
      - Draw number of events from Poisson(lambda), where lambda = mean_p_rec * length
      - For each event:
          * Sample a valid GC site
          * Pick one chromatid from each sister pair at random
          * Apply gene conversion (donor -> initiator)
          * 50:50: crossover or not

    Modifies post_rep_cell in place. Returns it for convenience.
    """
    chroms_A = post_rep_cell["A"]
    chroms_B = post_rep_cell["B"]

    # Compute lambda from the A homolog (both share same p_rec)
    ref_chrom = chroms_A[0]
    mean_p_rec = sum((e - s + 1) * v for s, e, v in ref_chrom["p_rec"]) / ref_chrom["length"]
    lam = mean_p_rec * ref_chrom["length"]

    n_events = _poisson_draw(lam)

    for _ in range(n_events):
        # Pick one chromatid from each homolog pair
        initiator = random.choice(chroms_A)
        donor = random.choice(chroms_B)

        site = _sample_gc_site(initiator, gc_min, gc_max)
        if site is None:
            continue
        gc_start, gc_end = site

        _apply_gene_conversion(initiator, donor, gc_start, gc_end)

        if random.random() < 0.5:
            _apply_crossover(initiator, donor, gc_start, gc_end)

    return post_rep_cell


def _poisson_draw(lam):
    """Draw from Poisson(lam). Uses Knuth algorithm for small lam."""
    if lam <= 0:
        return 0
    L = math.exp(-lam)
    k, p = 0, 1.0
    while p > L:
        k += 1
        p *= random.random()
    return k - 1


# ---------------------------------------------------------------------------
# Segregation
# ---------------------------------------------------------------------------

def segregate(post_rep_cell):
    """
    Segregate a post-replication cell into two daughter pre-replication cells.

    For each homolog pair, one chromatid goes to daughter 1 and one to daughter 2.
    Selection is random (centromere-linked, independent per homolog pair).

    Returns: (daughter1, daughter2)
        Each daughter is {"A": chr, "B": chr}
    """
    def split_pair(pair):
        idx = random.randint(0, 1)
        return pair[idx], pair[1 - idx]

    a1, a2 = split_pair(post_rep_cell["A"])
    b1, b2 = split_pair(post_rep_cell["B"])

    daughter1 = {"A": a1, "B": b1}
    daughter2 = {"A": a2, "B": b2}
    return daughter1, daughter2


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------

def select(daughter1, daughter2):
    """Randomly select one daughter cell for the next generation."""
    return random.choice([daughter1, daughter2])


# ---------------------------------------------------------------------------
# Full generation cycle
# ---------------------------------------------------------------------------

def run_generation(pre_rep_cell, gc_min=100, gc_max=5000):
    """
    Run one complete cell cycle:
        replication -> recombination -> segregation -> selection

    Returns the selected daughter as a pre-replication cell.
    """
    post_rep = replicate_cell(pre_rep_cell)
    recombine(post_rep, gc_min=gc_min, gc_max=gc_max)
    d1, d2 = segregate(post_rep)
    return select(d1, d2)


def run_simulation(genome, n_gen=10, gc_min=100, gc_max=5000):
    """
    Run the simulation for n_gen generations starting from genome.

    Returns the final pre-replication cell.
    """
    cell = genome
    for gen in range(1, n_gen + 1):
        cell = run_generation(cell, gc_min=gc_min, gc_max=gc_max)
    return cell


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

HAPLOTYPE_COLORS = {
    "A": "#2A36D5",    # blue
    "B": "#E1341E",    # red
    "het": "#DDDDDD",  # light grey — heterozygous
}


def compute_loh(cell):
    """
    Compute LOH state across the chromosome by intersecting the two homologs.

    For each position:
        homolog A == homolog B  ->  that haplotype label (LOH)
        homolog A != homolog B  ->  "het"

    Returns a list of (start, end, state) tuples where state is "A", "B", or "het".
    """
    segs_a = cell["A"]["haplotype"]
    segs_b = cell["B"]["haplotype"]
    result = []

    # Walk both segment lists simultaneously using a two-pointer merge
    i, j = 0, 0
    while i < len(segs_a) and j < len(segs_b):
        a_start, a_end, a_val = segs_a[i]
        b_start, b_end, b_val = segs_b[j]

        # Intersection of current segments
        lo = max(a_start, b_start)
        hi = min(a_end, b_end)

        if lo <= hi:
            state = a_val if a_val == b_val else "het"
            result.append((lo, hi, state))

        # Advance whichever segment ends first
        if a_end <= b_end:
            i += 1
        else:
            j += 1

    return merge_segments(result)


def _draw_chromosome_bar(ax, segments, y, bar_height, color_map, chrom):
    """Draw one horizontal bar of colored segments plus a centromere tick."""
    for seg_start, seg_end, val in segments:
        color = color_map.get(val, "#999999")
        ax.barh(
            y, seg_end - seg_start + 1,
            left=seg_start - 1,
            height=bar_height,
            color=color,
            edgecolor="none"
        )
    cen_start, cen_end = _cen_interval(chrom)
    cen_mid = (cen_start + cen_end) / 2
    ax.plot([cen_mid, cen_mid], [y - bar_height / 2, y + bar_height / 2],
            color="black", linewidth=2, zorder=5)


def plot_cell(cell, title="Haplotype map", ax=None):
    """
    Draw a linear haplotype map of the current cell.

    cell: pre-replication dict {"A": chr, "B": chr}
    Each homolog is drawn as a horizontal bar; haplotype segments are colored.
    The centromere is marked with a vertical tick.
    """
    show = ax is None
    if show:
        fig, ax = plt.subplots(figsize=(10, 2))

    homologs = [("A", cell["A"]), ("B", cell["B"])]
    bar_height = 0.5

    for idx, (label, chrom) in enumerate(homologs):
        _draw_chromosome_bar(ax, chrom["haplotype"], idx, bar_height,
                             HAPLOTYPE_COLORS, chrom)

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["A", "B"])
    ax.set_xlabel("Position (bp)")
    ax.set_title(title)
    ax.set_xlim(0, cell["A"]["length"])

    legend_patches = [
        mpatches.Patch(color=HAPLOTYPE_COLORS["A"], label="Haplotype A"),
        mpatches.Patch(color=HAPLOTYPE_COLORS["B"], label="Haplotype B"),
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.15, 1), loc="upper right", fontsize=8)

    if show:
        plt.tight_layout()
        plt.show()


def plot_loh(cell, title="LOH map", ax=None):
    """
    Draw a single LOH track for the cell.

    Blue   = homozygous A (LOH-A)
    Orange = homozygous B (LOH-B)
    Grey   = heterozygous (no LOH)

    The centromere position is taken from homolog A (both share the same CEN).
    """
    show = ax is None
    if show:
        fig, ax = plt.subplots(figsize=(10, 1.2))

    loh_segs = compute_loh(cell)
    _draw_chromosome_bar(ax, loh_segs, 0, 0.5, HAPLOTYPE_COLORS, cell["A"])

    ax.set_yticks([0])
    ax.set_yticklabels(["LOH"])
    ax.set_xlabel("Position (bp)")
    ax.set_title(title)
    ax.set_xlim(0, cell["A"]["length"])
    ax.set_ylim(-0.5, 0.5)

    legend_patches = [
        mpatches.Patch(color=HAPLOTYPE_COLORS["A"],   label="Homozygous A"),
        mpatches.Patch(color=HAPLOTYPE_COLORS["B"],   label="Homozygous B"),
        mpatches.Patch(color=HAPLOTYPE_COLORS["het"], label="Heterozygous"),
    ]
    ax.legend(handles=legend_patches, bbox_to_anchor=(1.15, 1), loc="upper right", fontsize=8)

    if show:
        plt.tight_layout()
        plt.show()


def plot_cell_and_loh(cell, title="", savepath=None):
    """
    Combined figure: haplotype map (top) and LOH track (bottom).
    Pass savepath to write to a file instead of displaying interactively.
    """
    fig, axes = plt.subplots(
        2, 1,
        figsize=(12, 3.5),
        gridspec_kw={"height_ratios": [2, 1]},
        sharex=True
    )
    fig.subplots_adjust(hspace=0.35)

    hap_title = f"{title} — Haplotypes" if title else "Haplotypes"
    loh_title = f"{title} — LOH" if title else "LOH"

    plot_cell(cell, title=hap_title, ax=axes[0])
    plot_loh(cell,  title=loh_title,  ax=axes[1])

    axes[1].set_xlabel("Position (bp)")

    if savepath:
        plt.savefig(savepath, dpi=150, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Mitotic recombination simulator")
    parser.add_argument("genome_file", help="CSV with columns: chromosome, length, CEN")
    parser.add_argument("--n_gen", type=int, default=10, help="Number of generations")
    parser.add_argument("--p_rec", type=float, default=1e-5, help="Base recombination probability per bp")
    parser.add_argument("--gc_min", type=int, default=100, help="Min GC tract size (bp)")
    parser.add_argument("--gc_max", type=int, default=5000, help="Max GC tract size (bp)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--plot", action="store_true", help="Plot final haplotype and LOH state")
    parser.add_argument("--plot-out", type=str, default=None,
                        help="Save plot to this path (e.g. out.png) instead of displaying interactively")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    genome = load_genome(args.genome_file, p_rec_default=args.p_rec)
    print(f"Loaded chromosome '{genome['A']['name']}', length {genome['A']['length']} bp")
    print(f"Running {args.n_gen} generation(s)...")

    final_cell = run_simulation(genome, n_gen=args.n_gen, gc_min=args.gc_min, gc_max=args.gc_max)

    print("\nFinal haplotype state:")
    for homolog in ("A", "B"):
        print(f"  Homolog {homolog}:")
        for seg in final_cell[homolog]["haplotype"]:
            print(f"    {seg[0]:>8} – {seg[1]:>8}  [{seg[2]}]")

    print("\nLOH state:")
    for seg in compute_loh(final_cell):
        print(f"    {seg[0]:>8} – {seg[1]:>8}  [{seg[2]}]")

    if args.plot or args.plot_out:
        plot_cell_and_loh(
            final_cell,
            title=f"After {args.n_gen} generation(s)",
            savepath=args.plot_out
        )
