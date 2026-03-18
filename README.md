# recombination_sim

A Python simulation of mitotic recombination and loss of heterozygosity (LOH) in a diploid yeast chromosome. Models a complete cell cycle — replication, inter-homolog recombination, segregation, and selection — and produces linear haplotype and LOH maps of the resulting genome state.

---

## Features

- Diploid chromosome representation using compact run-length encoded interval tuples
- Inter-homolog gene conversion with configurable tract size
- 50:50 crossover / non-crossover outcome per recombination event
- Poisson-distributed recombination events per generation
- Centromere-protected interval (GC placement disallowed at CEN)
- GC tracts clipped to chromosome ends
- Centromere-linked segregation mimicking random assortment
- Combined haplotype + LOH visualization
- Reproducible runs via random seed
- Single dependency beyond the standard library: `matplotlib`

---

## Requirements

- Python ≥ 3.9
- `matplotlib`

```bash
pip install matplotlib
```

---

## Installation

No installation required. Clone or download the repository and run directly:

```bash
git clone https://github.com/yourlab/recombination_sim.git
cd recombination_sim
```

---

## Input format

A CSV file with three columns specifying the chromosome name, length in bp, and centromere midpoint position:

```
chromosome, length, CEN
I, 500000, 150200
```

In single-chromosome mode (current default), only the first row is used. An example input file is provided as `genome_chrI.csv`.

---

## Usage

### Command line

```bash
python simRec.py genome_chrI.csv [options]
```

| Argument | Type | Default | Description |
|---|---|---|---|
| `genome_file` | positional | — | Path to genome CSV |
| `--n_gen` | int | `10` | Number of generations to simulate |
| `--p_rec` | float | `1e-5` | Base recombination probability per bp per generation |
| `--gc_min` | int | `100` | Minimum gene conversion tract size (bp) |
| `--gc_max` | int | `5000` | Maximum gene conversion tract size (bp) |
| `--seed` | int | `None` | Random seed for reproducibility |
| `--plot` | flag | off | Display combined haplotype + LOH plot interactively |
| `--plot-out` | str | `None` | Save plot to file (PNG or PDF) instead of displaying |

### Examples

Run 20 generations at elevated recombination rate with a fixed seed:

```bash
python simRec.py genome_chrI.csv --n_gen 20 --p_rec 1e-4 --seed 42
```

Run and save a plot:

```bash
python simRec.py genome_chrI.csv --n_gen 20 --p_rec 1e-4 --seed 42 --plot-out results.png
```

### Example output

```
Loaded chromosome 'I', length 500000 bp
Running 20 generation(s)...

Final haplotype state:
  Homolog A:
         1 –     4285  [B]
      4286 –     8059  [A]
      8060 –     8630  [B]
         ...

LOH state:
         1 –     4285  [B]
      4286 –     8059  [A]
      8060 –     8630  [het]
         ...
```

---

## Visualization

The `--plot` / `--plot-out` flags produce a two-panel figure:

- **Top panel — Haplotype map:** one bar per homolog (A and B), colored by parental identity at each position. The centromere is marked with a vertical black tick.
- **Bottom panel — LOH map:** a single bar showing genomic intervals that are homozygous (LOH-A in blue, LOH-B in orange) or heterozygous (light grey).

![Example plot](example_plot.png)

---

## Python API

The simulation can be driven programmatically without the CLI:

```python
from simRec import load_genome, run_simulation, plot_cell_and_loh, compute_loh

# Load genome from file
genome = load_genome("genome_chrI.csv", p_rec_default=1e-4)

# Run simulation
final_cell = run_simulation(genome, n_gen=20, gc_min=100, gc_max=5000)

# Inspect LOH state
for start, end, state in compute_loh(final_cell):
    print(f"{start:>8} – {end:>8}  [{state}]")

# Plot
plot_cell_and_loh(final_cell, title="20 generations", savepath="out.png")
```

Individual cycle steps are also importable for experimentation:

```python
from simRec import (
    replicate_cell,   # S phase: duplicate homologs into sister chromatid pairs
    recombine,        # recombination: GC ± crossover on post-replication cell
    segregate,        # segregation: split into two daughter cells
    select,           # selection: randomly choose one daughter
)
```

---

## Data structures

### Chromosome object

A chromosome is a `dict` with the following slots. All interval data are stored as lists of `(start, end, value)` tuples (1-based, inclusive coordinates):

| Slot | Type | Description |
|---|---|---|
| `name` | `str` | Chromosome name (e.g. `"I"`) |
| `length` | `int` | Chromosome length in bp |
| `position` | `[(int, int, None)]` | Fixed span; scalar metadata |
| `cen` | `[(int, int, None)]` | Centromere interval; protected from GC placement |
| `p_rec` | `[(int, int, float)]` | Recombination probability per bp; supports sub-intervals for hotspots |
| `haplotype` | `[(int, int, str)]` | Parental identity per interval; values `"A"` or `"B"` |

After a gene conversion without crossover, a haplotype slot might look like:

```python
[(1, 100000, "A"), (100001, 100100, "B"), (100101, 500000, "A")]
```

### Cell objects

**Pre-replication / post-segregation:**
```python
cell = {"A": chr_A, "B": chr_B}
```

**Post-replication (during recombination):**
```python
cell = {
    "A": [chrA1, chrA2],   # sister chromatid pair from homolog A
    "B": [chrB1, chrB2],   # sister chromatid pair from homolog B
}
```

---

## Recombination model

Each generation proceeds as follows:

1. **Replication** — each homolog is duplicated into a sister chromatid pair, producing four chromatids total.
2. **Recombination** — the number of events is drawn from a Poisson distribution with λ = `p_rec` × chromosome length. For each event:
   - A gene conversion (GC) site is sampled from non-CEN sequence, weighted by the `p_rec` slot.
   - GC tract size is drawn uniformly from [`gc_min`, `gc_max`] bp, centered on the site, clipped to chromosome ends, and rejected if it overlaps the centromere.
   - One chromatid is chosen at random from each sister pair (inter-homolog only).
   - The GC tract is copied from the donor chromatid into the initiating chromatid.
   - With 50% probability, a crossover occurs: the telomere-distal arm relative to the GC tract is swapped between the two chromatids.
3. **Segregation** — one chromatid from each sister pair is assigned to each daughter cell (centromere-linked, independent per homolog pair).
4. **Selection** — one daughter is chosen at random to seed the next generation.

---

## Testing

```bash
pip install pytest
pytest test_simRec.py -v
```

30 tests covering interval helpers, GC placement and CEN exclusion, gene conversion, crossover directionality, replication independence, segregation correctness, LOH computation, and end-to-end simulation integrity.

---

## Planned extensions

- Multi-chromosome support (all 16 *S. cerevisiae* chromosomes)
- Recombination hotspots via sub-interval `p_rec` segments
- Sister chromatid exchange (SCE) with configurable `P_sister`
- Generational history tracking for LOH accumulation plots
- PyRanges backend option for large-scale performance

---

## License

MIT
