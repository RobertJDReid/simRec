"""
recombination_sim/test_simRec.py

Unit tests covering interval helpers, GC placement, crossover logic,
replication, segregation, and a smoke test for full simulation.
"""

import random
import pytest
from simRec import (
    make_chromosome, clone_chromosome,
    query_value, replace_interval, merge_segments,
    load_genome, replicate_cell, recombine, segregate, select,
    run_generation, run_simulation,
    _sample_gc_site, _apply_gene_conversion, _apply_crossover,
    _cen_interval, _poisson_draw,
)


# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

def test_query_value_basic():
    segs = [(1, 100, "A"), (101, 200, "B")]
    assert query_value(segs, 1) == "A"
    assert query_value(segs, 100) == "A"
    assert query_value(segs, 101) == "B"
    assert query_value(segs, 200) == "B"

def test_query_value_missing():
    with pytest.raises(ValueError):
        query_value([(1, 100, "A")], 101)

def test_replace_interval_full():
    segs = [(1, 500000, "A")]
    result = replace_interval(segs, 100001, 100100, "B")
    assert result == [(1, 100000, "A"), (100001, 100100, "B"), (100101, 500000, "A")]

def test_replace_interval_at_start():
    segs = [(1, 500000, "A")]
    result = replace_interval(segs, 1, 1000, "B")
    assert result[0] == (1, 1000, "B")
    assert result[1] == (1001, 500000, "A")

def test_replace_interval_at_end():
    segs = [(1, 500000, "A")]
    result = replace_interval(segs, 499000, 500000, "B")
    assert result[-1] == (499000, 500000, "B")

def test_merge_segments_adjacent_same():
    segs = [(1, 100, "A"), (101, 200, "A"), (201, 300, "B")]
    assert merge_segments(segs) == [(1, 200, "A"), (201, 300, "B")]

def test_merge_segments_no_merge():
    segs = [(1, 100, "A"), (101, 200, "B")]
    assert merge_segments(segs) == segs

def test_replace_then_same_value_merges():
    # Replacing a middle segment with the same value as neighbors should merge
    segs = [(1, 500000, "A")]
    result = replace_interval(segs, 100, 200, "A")
    assert result == [(1, 500000, "A")]


# ---------------------------------------------------------------------------
# Chromosome construction
# ---------------------------------------------------------------------------

def test_make_chromosome_slots():
    chrom = make_chromosome("I", 500000, 150100, 150300, p_rec_default=1e-5, parent="A")
    assert chrom["length"] == 500000
    assert chrom["haplotype"] == [(1, 500000, "A")]
    assert chrom["cen"] == [(150100, 150300, None)]
    assert chrom["p_rec"] == [(1, 500000, 1e-5)]

def test_clone_independence():
    chrom = make_chromosome("I", 500000, 150100, 150300)
    clone = clone_chromosome(chrom)
    clone["haplotype"] = replace_interval(clone["haplotype"], 1, 1000, "B")
    assert chrom["haplotype"] == [(1, 500000, "A")]


# ---------------------------------------------------------------------------
# GC site sampling
# ---------------------------------------------------------------------------

def _simple_chrom(parent="A"):
    return make_chromosome("I", 500000, 150100, 150300, p_rec_default=1e-5, parent=parent)

def test_gc_site_within_bounds():
    random.seed(42)
    chrom = _simple_chrom()
    none_count = 0
    for _ in range(200):
        site = _sample_gc_site(chrom, gc_min=100, gc_max=5000)
        if site is None:
            none_count += 1
            continue
        gc_start, gc_end = site
        assert gc_start >= 1
        assert gc_end <= chrom["length"]
        assert gc_start <= gc_end
    # On a 500 kb chromosome with a tiny CEN, None should essentially never occur
    assert none_count == 0, f"Unexpected None placements: {none_count}"

def test_gc_site_no_cen_overlap():
    random.seed(7)
    chrom = _simple_chrom()
    cen_start, cen_end = _cen_interval(chrom)
    for _ in range(100):
        site = _sample_gc_site(chrom, gc_min=100, gc_max=5000)
        if site is None:
            continue
        gc_start, gc_end = site
        overlap = gc_end >= cen_start and gc_start <= cen_end
        assert not overlap, f"GC {gc_start}-{gc_end} overlaps CEN {cen_start}-{cen_end}"


# ---------------------------------------------------------------------------
# Gene conversion
# ---------------------------------------------------------------------------

def test_gene_conversion_copies_donor():
    initiator = _simple_chrom("A")
    donor = _simple_chrom("B")
    _apply_gene_conversion(initiator, donor, 100001, 100100)
    assert query_value(initiator["haplotype"], 100050) == "B"
    assert query_value(initiator["haplotype"], 1) == "A"
    assert query_value(initiator["haplotype"], 500000) == "A"

def test_gene_conversion_does_not_touch_donor():
    initiator = _simple_chrom("A")
    donor = _simple_chrom("B")
    _apply_gene_conversion(initiator, donor, 100001, 100100)
    assert donor["haplotype"] == [(1, 500000, "B")]


# ---------------------------------------------------------------------------
# Crossover
# ---------------------------------------------------------------------------

def test_crossover_right_arm():
    """GC on right arm: region after GC should swap."""
    a = _simple_chrom("A")
    b = _simple_chrom("B")
    # GC on right arm (after centromere at 150100-150300)
    gc_start, gc_end = 200000, 201000
    _apply_gene_conversion(a, b, gc_start, gc_end)
    _apply_crossover(a, b, gc_start, gc_end)
    # After crossover, a should have "B" from gc_end+1 to end
    assert query_value(a["haplotype"], gc_end + 1) == "B"
    assert query_value(a["haplotype"], 500000) == "B"
    # And b should have "A" in that region
    assert query_value(b["haplotype"], gc_end + 1) == "A"

def test_crossover_left_arm():
    """GC on left arm (before centromere): region before GC should swap."""
    a = _simple_chrom("A")
    b = _simple_chrom("B")
    gc_start, gc_end = 10000, 10500
    _apply_gene_conversion(a, b, gc_start, gc_end)
    _apply_crossover(a, b, gc_start, gc_end)
    # Before the GC, a should now be "B"
    assert query_value(a["haplotype"], 1) == "B"
    assert query_value(a["haplotype"], gc_start - 1) == "B"
    assert query_value(b["haplotype"], 1) == "A"


# ---------------------------------------------------------------------------
# Replication
# ---------------------------------------------------------------------------

def test_replication_produces_four_chromatids():
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    post = replicate_cell(genome)
    assert len(post["A"]) == 2
    assert len(post["B"]) == 2

def test_replication_independence():
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    post = replicate_cell(genome)
    # Mutate one sister; the other should be unaffected
    post["A"][0]["haplotype"] = replace_interval(post["A"][0]["haplotype"], 1, 1000, "B")
    assert query_value(post["A"][1]["haplotype"], 500) == "A"


# ---------------------------------------------------------------------------
# Segregation
# ---------------------------------------------------------------------------

def test_segregation_produces_two_daughters():
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    post = replicate_cell(genome)
    d1, d2 = segregate(post)
    assert "A" in d1 and "B" in d1
    assert "A" in d2 and "B" in d2

def test_segregation_each_chromatid_in_exactly_one_daughter():
    random.seed(0)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    post = replicate_cell(genome)
    chrA1_id = id(post["A"][0])
    chrA2_id = id(post["A"][1])
    d1, d2 = segregate(post)
    d1A_id = id(d1["A"])
    d2A_id = id(d2["A"])
    assert {d1A_id, d2A_id} == {chrA1_id, chrA2_id}


# ---------------------------------------------------------------------------
# Poisson draw
# ---------------------------------------------------------------------------

def test_poisson_draw_zero_lambda():
    assert _poisson_draw(0) == 0

def test_poisson_draw_reasonable_mean():
    random.seed(1)
    samples = [_poisson_draw(3.0) for _ in range(2000)]
    mean = sum(samples) / len(samples)
    assert abs(mean - 3.0) < 0.2


# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------

def test_full_simulation_runs():
    random.seed(99)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    final = run_simulation(genome, n_gen=5)
    assert "A" in final and "B" in final
    # Haplotype segments should cover the full chromosome without gaps
    for homolog in ("A", "B"):
        segs = final[homolog]["haplotype"]
        assert segs[0][0] == 1
        assert segs[-1][1] == 500000
        for i in range(len(segs) - 1):
            assert segs[i][1] + 1 == segs[i + 1][0], "Gap in haplotype segments"

def test_haplotype_values_valid():
    random.seed(42)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    final = run_simulation(genome, n_gen=20)
    for homolog in ("A", "B"):
        for _, _, val in final[homolog]["haplotype"]:
            assert val in ("A", "B")


# ---------------------------------------------------------------------------
# LOH computation
# ---------------------------------------------------------------------------

from simRec import compute_loh

def test_loh_fully_heterozygous():
    """Unrecombined cell: A homolog is all-A, B homolog is all-B -> all het."""
    cell = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    loh = compute_loh(cell)
    assert loh == [(1, 500000, "het")]

def test_loh_fully_homozygous_A():
    """Both homologs all-A -> all LOH-A."""
    cell = {"A": _simple_chrom("A"), "B": _simple_chrom("A")}
    loh = compute_loh(cell)
    assert loh == [(1, 500000, "A")]

def test_loh_fully_homozygous_B():
    """Both homologs all-B -> all LOH-B."""
    cell = {"A": _simple_chrom("B"), "B": _simple_chrom("B")}
    loh = compute_loh(cell)
    assert loh == [(1, 500000, "B")]

def test_loh_partial_conversion():
    """
    Homolog A has a B-converted region 100001-100100.
    Homolog B is all-B.
    Expect: 1-100000 het, 100001-100100 LOH-B, 100101-500000 het.
    """
    cell = {
        "A": _simple_chrom("A"),
        "B": _simple_chrom("B"),
    }
    cell["A"]["haplotype"] = replace_interval(cell["A"]["haplotype"], 100001, 100100, "B")
    loh = compute_loh(cell)
    assert (100001, 100100, "B") in loh
    # flanking regions should be het
    het_segs = [(s, e, v) for s, e, v in loh if v == "het"]
    assert any(s <= 1 and e >= 100000 for s, e, _ in het_segs)
    assert any(s <= 100101 and e >= 500000 for s, e, _ in het_segs)

def test_loh_covers_full_chromosome():
    """LOH segments should cover 1 to length with no gaps."""
    random.seed(5)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    final = run_simulation(genome, n_gen=10)
    loh = compute_loh(final)
    assert loh[0][0] == 1
    assert loh[-1][1] == 500000
    for i in range(len(loh) - 1):
        assert loh[i][1] + 1 == loh[i + 1][0], "Gap in LOH segments"

def test_loh_values_valid():
    random.seed(13)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    final = run_simulation(genome, n_gen=15)
    for _, _, val in compute_loh(final):
        assert val in ("A", "B", "het")
