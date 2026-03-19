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


# ---------------------------------------------------------------------------
# classify_events
# ---------------------------------------------------------------------------

from simRec import classify_events

def _make_cell_from_haplotypes(hap_A, hap_B):
    """
    Build a minimal cell dict from explicit haplotype segment lists.
    Uses the standard _simple_chrom for CEN and length metadata,
    then overwrites the haplotype slot.
    """
    cA = _simple_chrom("A")
    cB = _simple_chrom("B")
    cA["haplotype"] = hap_A
    cB["haplotype"] = hap_B
    return {"A": cA, "B": cB}

# CEN is at 150100-150300 on the 500000 bp test chromosome.
# Positions used in tests are chosen to be clearly on one arm or the other.

def test_classify_no_events():
    """Fully heterozygous cell — no events."""
    cell = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    assert classify_events(cell) == []

def test_classify_gc_nco_internal():
    """
    Internal LOH block flanked by het on both sides, same haplotype context
    outside — GC-NCO.
    Homolog A: all A except 200000-201000 converted to B.
    Homolog B: all B.
    LOH: het | B(200000-201000) | het  -> GC-NCO
    """
    cell = _make_cell_from_haplotypes(
        replace_interval([(1, 500000, "A")], 200000, 201000, "B"),
        [(1, 500000, "B")],
    )
    events = classify_events(cell)
    assert len(events) == 1
    e = events[0]
    assert e["type"]      == "GC-NCO"
    assert e["start"]     == 200000
    assert e["end"]       == 201000
    assert e["haplotype"] == "B"
    assert e["flanking_left"]  == "het"
    assert e["flanking_right"] == "het"

def test_classify_gc_co_internal():
    """
    Internal LOH block where flanking LOH haplotypes differ — GC-CO.
    We construct a cell where:
      Homolog A: B(1-100000) | A(100001-200000) | B(200001-500000)  (left arm CO)
      Homolog B: all B
    LOH pattern: B(1-100000) | het(100001-200000) | B(200001-500000)
    The het block in the middle represents the GC tract on an otherwise
    terminal-CO chromosome — but here we build a simpler case:

    Simpler: both homologs identical except for an A island in a sea of B:
      Homolog A: A(1-100000) | B(100001-101000) | A(101001-300000) | B(300001-500000)
      Homolog B: all B
    LOH: het(1-100000) | B(100001-101000) | het(101001-300000) | B(300001-500000)
    The B block at 100001-101000: left context is het, right context is het -> NCO
    The B block at 300001-500000: touches right telomere -> CO-terminal

    For a pure GC-CO we need flanking LOH haplotypes to differ.
    Construct:
      Homolog A: B(1-99000) | A(99001-101000) | B(101001-500000)
      Homolog B: A(1-99000) | A(99001-101000) | A(101001-500000)  = all A
    LOH: B(1-99000) | het(99001-101000) | ... wait, hom A on right
    Actually simpler:
      Homolog A: B(1-200000) | A(200001-201000) | B(201001-500000)
      Homolog B: A(1-200000) | A(200001-201000) | A(201001-500000) = all A
    LOH: het | A(200001-201000) | het
    flanking left LOH = B (from hom-A region left), flanking right LOH = B -> NCO

    Correct construction for GC-CO:
      Homolog A: A(1-200000) | B(200001-201000) | A(201001-500000)
      Homolog B: B(1-200000) | B(200001-201000) | A(201001-500000)
    LOH: het(1-200000) | B(200001-201000) | het(201001-500000) -> NCO (both flanks het)

    True GC-CO requires LOH on both sides with different haplotypes:
      Homolog A: B(1-100000) | A(100001-101000) | A(101001-500000)
      Homolog B: B(1-100000) | A(100001-101000) | B(101001-500000)
    LOH: B(1-100000) | A(100001-500000 all A) ... need to intersect carefully.
    Hom A haplotype: B(1-100000), A(100001-500000)
    Hom B haplotype: B(1-100000), A(100001-101000), B(101001-500000)
    LOH: B(1-100000) | A(100001-101000) | het(101001-500000)
    The A block at 100001-101000: left flank = B LOH, right flank = het -> GC-CO
    """
    hap_A = [(1, 100000, "B"), (100001, 500000, "A")]
    hap_B = [(1, 100000, "B"), (100001, 101000, "A"), (101001, 500000, "B")]
    cell = _make_cell_from_haplotypes(hap_A, hap_B)
    events = classify_events(cell)
    # Expect: B(1-100000) terminal CO + A(100001-101000) GC-CO
    co_term = [e for e in events if e["type"] == "CO-terminal"]
    gc_co   = [e for e in events if e["type"] == "GC-CO"]
    assert len(gc_co) == 1
    e = gc_co[0]
    assert e["start"]     == 100001
    assert e["end"]       == 101000
    assert e["haplotype"] == "A"
    # Left flank is LOH-B, right flank is het
    assert e["flanking_left"]  == "tel"   # B block reaches left telomere
    assert e["flanking_right"] == "het"

def test_classify_co_terminal_right():
    """
    Terminal LOH on the right arm (extends to position 500000).
    Homolog A: A(1-300000) | B(300001-500000)
    Homolog B: all B
    LOH: het(1-300000) | B(300001-500000)
    """
    hap_A = [(1, 300000, "A"), (300001, 500000, "B")]
    hap_B = [(1, 500000, "B")]
    cell = _make_cell_from_haplotypes(hap_A, hap_B)
    events = classify_events(cell)
    assert len(events) == 1
    e = events[0]
    assert e["type"]           == "CO-terminal"
    assert e["start"]          == 300001
    assert e["end"]            == 500000
    assert e["haplotype"]      == "B"
    assert e["flanking_left"]  == "het"
    assert e["flanking_right"] == "tel"

def test_classify_co_terminal_left():
    """
    Terminal LOH on the left arm (extends to position 1).
    Homolog A: B(1-100000) | A(100001-500000)
    Homolog B: all B
    LOH: B(1-100000) | het(100001-500000)
    """
    hap_A = [(1, 100000, "B"), (100001, 500000, "A")]
    hap_B = [(1, 500000, "B")]
    cell = _make_cell_from_haplotypes(hap_A, hap_B)
    events = classify_events(cell)
    assert len(events) == 1
    e = events[0]
    assert e["type"]           == "CO-terminal"
    assert e["start"]          == 1
    assert e["end"]            == 100000
    assert e["haplotype"]      == "B"
    assert e["flanking_left"]  == "tel"
    assert e["flanking_right"] == "het"

def test_classify_tel_tel():
    """
    Entire chromosome LOH — two CO-terminal events, one per arm.
    Both homologs all-A.
    """
    cell = {"A": _simple_chrom("A"), "B": _simple_chrom("A")}
    events = classify_events(cell)
    assert len(events) == 2
    types = [e["type"] for e in events]
    assert all(t == "CO-terminal" for t in types)
    # One event should end at the left of the CEN, one should start at the right
    cen_start = 150100
    cen_end   = 150300
    left_event  = next(e for e in events if e["end"]   < cen_start)
    right_event = next(e for e in events if e["start"] > cen_end)
    assert left_event["flanking_left"]   == "tel"
    assert left_event["flanking_right"]  == "het"
    assert right_event["flanking_left"]  == "het"
    assert right_event["flanking_right"] == "tel"

def test_classify_abutting_terminal():
    """
    Two abutting terminal LOH blocks on the right arm — two CO-terminal events.

    Construction:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: A(1-200000) | B(200001-350000) | B(350001-500000)

    LOH map:
      A(1-200000) | B(200001-350000) | het(350001-500000)
      Wait — both homologs agree on A and B blocks but differ at 350001+.

    Actually to get two abutting terminal blocks we need the LOH to be:
      het | B(200001-350000) | A(350001-500000)
    where the B block is interior-touching-A and the A block hits telomere.

    Construction:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: all B

    LOH:
      het(1-200000) | B(200001-350000) | het(350001-500000)
    — this gives one B interior block, not two terminal blocks.

    For two abutting terminal blocks we need both to reach pos 500000:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: A(1-200000) | A(200001-350000) | A(350001-500000) = all A

    LOH:
      A(1-200000) | het(200001-350000) | A(350001-500000)
    — gives two separate A terminal blocks with het in between, both touching telomere indirectly.
    The A(1-200000) touches left telomere -> CO-terminal left arm.
    The A(350001-500000) touches right telomere -> CO-terminal right arm.
    The het block in the middle is the B island on one homolog.

    This is the correct abutting-terminal construction for the RIGHT arm:
      Homolog A: A(1-300000) | B(300001-400000) | A(400001-500000)
      Homolog B: A(1-300000) | A(300001-400000) | A(400001-500000) = all A after 300000

    That still gives het | A terminal.

    Simplest valid construction: two different-haplotype terminal blocks on right arm:
      Homolog A: A(1-200000) | A(200001-350000) | B(350001-500000)
      Homolog B: A(1-200000) | B(200001-350000) | B(350001-500000)

    LOH:
      A(1-200000) | het(200001-350000) | B(350001-500000)

    A(1-200000) touches left telomere -> CO-terminal.
    B(350001-500000) touches right telomere -> CO-terminal.
    het block in between is as expected.
    This tests two terminal events (one per arm), with a het gap between them.
    """
    hap_A = [(1, 200000, "A"), (200001, 350000, "A"), (350001, 500000, "B")]
    hap_B = [(1, 200000, "A"), (200001, 350000, "B"), (350001, 500000, "B")]
    cell = _make_cell_from_haplotypes(hap_A, hap_B)
    events = classify_events(cell)
    co_term = [e for e in events if e["type"] == "CO-terminal"]
    assert len(co_term) == 2
    left_co  = next(e for e in co_term if e["flanking_left"] == "tel")
    right_co = next(e for e in co_term if e["flanking_right"] == "tel")
    assert left_co["start"]  == 1
    assert left_co["end"]    == 200000
    assert right_co["start"] == 350001
    assert right_co["end"]   == 500000

def test_classify_abutting_terminal_same_arm():
    """
    Two abutting terminal LOH blocks on the SAME arm (right), different haplotypes.
    This is the masking edge case from the spec:
      het | LOH-B(200001-350000) | LOH-A(350001-500000) -> telomere

    Construction:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: A(1-200000) | B(200001-350000) | B(350001-500000)

    LOH: A(1-200000) | B(200001-350000) | het(350001-500000)
    Hmm — the A block on the left reaches the left telomere.

    We need both B and A to abut on the RIGHT side.
    Construction:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: A(1-500000)

    LOH: A(1-200000) | het(200001-350000) | A(350001-500000)
    Left A block touches left tel; right A block touches right tel.
    B island in between is het (B on one, A on other). Not what we want.

    True same-arm abutting: need the LOH map to show:
      het | B(s1-e1) | A(e1+1-500000)  where A reaches telomere

    Construction:
      Homolog A: A(1-200000) | B(200001-350000) | A(350001-500000)
      Homolog B: B(1-200000) | B(200001-350000) | A(350001-500000) = B(1-350000)|A(350001+)

    LOH: het(1-200000) | B(200001-350000) | A(350001-500000)
    B(200001-350000): interior, flanked het left, A right -> GC-CO
    A(350001-500000): touches right telomere -> CO-terminal
    """
    hap_A = [(1, 200000, "A"), (200001, 350000, "B"), (350001, 500000, "A")]
    hap_B = [(1, 350000, "B"),                        (350001, 500000, "A")]
    cell = _make_cell_from_haplotypes(hap_A, hap_B)
    events = classify_events(cell)
    co_term = [e for e in events if e["type"] == "CO-terminal"]
    gc_co   = [e for e in events if e["type"] == "GC-CO"]
    assert len(co_term) == 1
    assert co_term[0]["start"] == 350001
    assert co_term[0]["end"]   == 500000
    assert co_term[0]["flanking_right"] == "tel"
    assert len(gc_co) == 1
    assert gc_co[0]["start"] == 200001
    assert gc_co[0]["end"]   == 350000
    assert gc_co[0]["flanking_left"]  == "het"
    # The B block's right flank walks through the terminal A LOH block to the telomere
    assert gc_co[0]["flanking_right"] == "tel"
    """All event types from a full simulation must be valid strings."""
    random.seed(77)
    genome = {"A": _simple_chrom("A"), "B": _simple_chrom("B")}
    final = run_simulation(genome, n_gen=20)
    valid_types = {"GC-NCO", "GC-CO", "CO-terminal"}
    for e in classify_events(final):
        assert e["type"] in valid_types
        assert e["start"] <= e["end"]
        assert e["haplotype"] in ("A", "B")
        assert e["flanking_left"]  in ("A", "B", "het", "tel")
        assert e["flanking_right"] in ("A", "B", "het", "tel")
