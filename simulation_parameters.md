# simRec — Simulation Overview and Parameters

## Overview

simRec simulates mitotic recombination and loss of heterozygosity (LOH) in a diploid yeast chromosome across multiple generations. Each generation models a complete cell cycle: replication of both parental homologs into sister chromatid pairs, inter-homolog recombination via gene conversion with or without crossover, centromere-linked segregation into two daughter cells, and random selection of one daughter to seed the next generation. The resulting haplotype state is analyzed post-hoc using rules that mirror the interpretation of experimental LOH data, classifying observed LOH intervals as gene conversion non-crossover events (GC-NCO), gene conversion crossover events (GC-CO), or terminal crossovers (CO-terminal).

## Parameters

### Genome

- **Chromosome name** — identifier read from the input CSV (e.g. `I`)
- **Chromosome length** — total length in bp read from the input CSV
- **Centromere position** — midpoint of the centromere interval in bp, read from the input CSV
- **Centromere width** — size of the protected centromere interval in bp; gene conversion tracts are disallowed from overlapping this region (default: 200 bp, centered on the midpoint)

### Simulation

- **`--n-gen`** — number of generations to simulate (default: 10)
- **`--p-rec`** — base recombination probability per bp per generation; used to set a uniform `p_rec` value across the chromosome, which scales the Poisson-distributed number of recombination events per generation as λ = `p_rec` × chromosome length (default: 1 × 10⁻⁵)
- **`--gc-min`** — minimum gene conversion tract size in bp; tracts are drawn uniformly from [`gc-min`, `gc-max`] and centered on the recombination site (default: 100 bp)
- **`--gc-max`** — maximum gene conversion tract size in bp (default: 5000 bp)
- **`--seed`** — integer random seed for reproducibility; if not set, results will differ between runs

### Recombination model (fixed)

- **Partner selection** — inter-homolog only; sister chromatid exchange is not currently modeled
- **Crossover probability** — 50% per recombination event, independent of position
- **Crossover arm** — the telomere-distal segment relative to the gene conversion tract is swapped; which arm is determined by whether the tract falls left or right of the centromere
- **Centromere exclusion** — gene conversion tracts that would overlap the centromere interval are rejected and resampled
- **Tract clipping** — tracts extending beyond chromosome ends are clipped to the boundary

### Segregation (fixed)

- **Segregation mode** — centromere-linked random assortment; one chromatid from each sister pair is assigned independently to each daughter cell
- **Selection** — one of the two daughter cells is chosen at random to seed the next generation (neutral drift, no fitness model)
