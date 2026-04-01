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

def load_genome(filepath, p_rec_default=1e-5):
    with open(filepath) as fh:
        reader = csv.DictReader(fh, skipinitialspace=True)
        rows = list(reader)

    if not rows:
        raise ValueError("Genome file is empty.")

    row = rows[0]
    name      = row["chromosome"].strip()
    length    = int(row["length"])
    cen_start = int(row["CEN_start"])          # ← read directly
    cen_end   = int(row["CEN_end"])            # ← read directly

    if not (1 <= cen_start <= cen_end <= length):
        raise ValueError(
            f"Invalid CEN range [{cen_start}, {cen_end}] for chromosome "
            f"'{name}' of length {length}."
        )

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

def recombine(post_rep_cell, gc_min=100, gc_max=5000, co_prob=0.5):
    """
    Simulate recombination on a post-replication cell.

    For each chromosome:
      - Draw number of events from Poisson(lambda), where lambda = mean_p_rec * length
      - For each event:
          * Sample a valid GC site
          * Pick one chromatid from each sister pair at random
          * Apply gene conversion (donor -> initiator)
          * Crossover with probability co_prob, otherwise NCO

    co_prob : float in [0, 1]
        Probability that a recombination event results in a crossover.
        Default 0.5 (1:1 CO:NCO ratio).

    Modifies post_rep_cell in place.

    Returns
    -------
    list of dict
        One entry per recombination event that fired, with keys:
            chrom    : chromosome name
            type     : "GC" (gene conversion only) or "CO" (crossover + GC)
            start    : GC tract start position
            end      : GC tract end position
        Empty list if no events occurred.
    """
    chroms_A = post_rep_cell["A"]
    chroms_B = post_rep_cell["B"]

    # Compute lambda from the A homolog (both share same p_rec)
    ref_chrom = chroms_A[0]
    mean_p_rec = sum((e - s + 1) * v for s, e, v in ref_chrom["p_rec"]) / ref_chrom["length"]
    lam = mean_p_rec * ref_chrom["length"]

    n_events = _poisson_draw(lam)
    logged = []

    for _ in range(n_events):
        # Pick one chromatid from each homolog pair, randomly assigning
        # which homolog initiates and which donates
        pair1, pair2 = random.sample([chroms_A, chroms_B], 2)
        initiator = random.choice(pair1)
        donor     = random.choice(pair2)

        site = _sample_gc_site(initiator, gc_min, gc_max)
        if site is None:
            continue
        gc_start, gc_end = site

        converts_to = query_value(donor["haplotype"], (gc_start + gc_end) // 2)
        _apply_gene_conversion(initiator, donor, gc_start, gc_end)

        if random.random() < co_prob:
            _apply_crossover(initiator, donor, gc_start, gc_end)
            event_type = "CO"
        else:
            event_type = "NCO"

        logged.append({
            "chrom":       initiator["name"],
            "type":        event_type,
            "start":       gc_start,
            "end":         gc_end,
            "converts_to": converts_to,
        })

    return post_rep_cell, logged


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

def run_generation(pre_rep_cell, gc_min=100, gc_max=5000, co_prob=0.5):
    """
    Run one complete cell cycle:
        replication -> recombination -> segregation -> selection

    Returns
    -------
    cell : dict
        The selected daughter as a pre-replication cell.
    logged : list of dict
        Recombination events that fired this generation, each with keys:
        chrom, type ("GC" or "CO"), start, end.
    """
    post_rep = replicate_cell(pre_rep_cell)
    post_rep, logged = recombine(post_rep, gc_min=gc_min, gc_max=gc_max,
                                 co_prob=co_prob)
    d1, d2 = segregate(post_rep)
    return select(d1, d2), logged


def run_simulation(genome, n_gen=10, gc_min=100, gc_max=5000, co_prob=0.5):
    """
    Run the simulation for n_gen generations starting from genome.

    Returns
    -------
    cell : dict
        The final pre-replication cell.
    event_log : list of dict
        All recombination events across all generations, each with keys:
        gen (1-based), chrom, type ("GC" or "CO"), start, end.
    haplotype_snapshots : list of (int, dict)
        One (gen, cell) pair for each generation in which at least one
        recombination event fired, capturing the post-selection haplotype
        state. Generations with no events are omitted.
    """
    cell                = genome
    event_log           = []
    haplotype_snapshots = []
    for gen in range(1, n_gen + 1):
        cell, logged = run_generation(cell, gc_min=gc_min, gc_max=gc_max,
                                      co_prob=co_prob)
        if logged:
            for entry in logged:
                entry["gen"] = gen
                event_log.append(entry)
            haplotype_snapshots.append((gen, cell))
    return cell, event_log, haplotype_snapshots


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


def classify_events(cell, gc_max=5000):
    """
    Post-hoc classification of recombination events from the LOH pattern
    of a pre-replication cell.

    Classification rules
    --------------------
    The LOH segment list (from compute_loh) is scanned for non-"het" blocks.
    Each is classified as one of:

        NCO-GC       : Internal LOH block ≤ gc_max where the haplotype phase
                       on homolog A does NOT switch across the block.
 
        CO-GC        : Internal LOH block ≤ gc_max where the haplotype phase
                       on homolog A DOES switch across the block.
 
        DCO-internal : Internal LOH block > gc_max flanked by het on both sides.
                       Cannot be a gene conversion; interpreted as a double
                       crossover bracketing the region.
 
        CO-terminal  : LOH block that extends to a telomere (start == 1 or
                       end == chrom_length).
 
        TEL-TEL      : Entire chromosome is LOH. Reported as two CO-terminal
                       events anchored at the CEN boundaries.

    Blocks > gc_max with a LOH neighbor on either side receive type "NCO-GC"
    with complex=True, flagging them for case-by-case examination.

    Parameters
    ----------
    cell : dict
        Pre-replication cell {"A": chrom, "B": chrom}.
    gc_max : int
        Maximum gene conversion tract size used in the simulation (default 5000).
        Internal blocks exceeding this size are classified as DCO-internal or complex.

    Returns
    -------
    list of dict, one record per classified event, in chromosome order.
    """
    chrom        = cell["A"]
    chrom_length = chrom["length"]
    cen_start, cen_end = _cen_interval(chrom)

    hap_A    = cell["A"]["haplotype"]
    hap_B    = cell["B"]["haplotype"]
    loh_segs = compute_loh(cell)
    n        = len(loh_segs)
    events   = []

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _flank_context(idx, direction):
        """
        Return "het" or "tel" for the immediate context on one side of
        loh_segs[idx]. direction is "left" or "right".
        Only steps one segment — does NOT walk through adjacent LOH blocks.
        """
        if direction == "left":
            if loh_segs[idx][0] == 1:
                return "tel"
            neighbor_idx = idx - 1
        else:
            if loh_segs[idx][1] == chrom_length:
                return "tel"
            neighbor_idx = idx + 1

        if 0 <= neighbor_idx < n:
            nv = loh_segs[neighbor_idx][2]
            if nv == "het":
                return "het"
            # Neighbor is another LOH block — still return "het" as the
            # immediate separator context is the boundary itself
            return "het"
        return "tel"

    def _hap_A_at(pos):
        """Query homolog A haplotype at a given position."""
        return query_value(hap_A, pos)

    def _is_terminal_block(idx):
        """True if this LOH block touches a telomere."""
        s, e, v = loh_segs[idx]
        return v != "het" and (s == 1 or e == chrom_length)

    def _neighbor_loh_value(idx, direction):
        """
        Return the LOH state of the immediately neighboring segment.
        direction: "left" (idx-1) or "right" (idx+1).
        Returns "tel" if at chromosome boundary, otherwise the neighbor's state
        ("het", "A", or "B").
        """
        if direction == "left":
            if loh_segs[idx][0] == 1:
                return "tel"
            neighbor_idx = idx - 1
        else:
            if loh_segs[idx][1] == chrom_length:
                return "tel"
            neighbor_idx = idx + 1
        if 0 <= neighbor_idx < n:
            return loh_segs[neighbor_idx][2]
        return "tel"

    def _adjacent_to_terminal(fl, fr):
        """
        True if either immediate neighbor is a non-het LOH block (meaning this
        block is directly part of a continuous LOH cluster that reaches a telomere)
        OR if either flank is "tel" — but terminal blocks are handled separately,
        so here we only need to check for LOH-valued neighbors.
        A neighbor value of "A" or "B" means the next segment is itself LOH,
        placing this block inside a terminal LOH cluster.
        """
        return fl in ("A", "B") or fr in ("A", "B")

    # ------------------------------------------------------------------
    # Special case: TEL-TEL (no het anywhere)
    # ------------------------------------------------------------------
    loh_states = {v for _, _, v in loh_segs}
    if "het" not in loh_states:
        left_segs  = [(s, e, v) for s, e, v in loh_segs if e < cen_start]
        right_segs = [(s, e, v) for s, e, v in loh_segs if s > cen_end]
        if not left_segs:
            left_segs  = [(s, min(e, cen_start - 1), v)
                          for s, e, v in loh_segs if s < cen_start]
        if not right_segs:
            right_segs = [(max(s, cen_end + 1), e, v)
                          for s, e, v in loh_segs if e > cen_end]
        if left_segs:
            events.append({
                "type":                 "CO-terminal",
                "chrom":                chrom["name"],
                "start":                1,
                "end":                  cen_start - 1,
                "haplotype":            left_segs[-1][2],
                "flanking_left":        "tel",
                "flanking_right":       "het",
                "adjacent_to_terminal": False,
            })
        if right_segs:
            events.append({
                "type":                 "CO-terminal",
                "chrom":                chrom["name"],
                "start":                cen_end + 1,
                "end":                  chrom_length,
                "haplotype":            right_segs[0][2],
                "flanking_left":        "het",
                "flanking_right":       "tel",
                "adjacent_to_terminal": False,
            })
        return events

    # ------------------------------------------------------------------
    # General case
    # ------------------------------------------------------------------
    for i in range(n):
        seg_start, seg_end, state = loh_segs[i]

        if state == "het":
            continue

        hap = state
        touches_left_tel  = (seg_start == 1)
        touches_right_tel = (seg_end   == chrom_length)

        # Flanking context: actual LOH state of the immediate neighbor segment,
        # or "tel" if at a chromosome boundary.
        fl = _neighbor_loh_value(i, "left")
        fr = _neighbor_loh_value(i, "right")

        if touches_left_tel or touches_right_tel:
            events.append({
                "type":                 "CO-terminal",
                "chrom":                chrom["name"],
                "start":                seg_start,
                "end":                  seg_end,
                "haplotype":            hap,
                "flanking_left":        fl,
                "flanking_right":       fr,
                "adjacent_to_terminal": False,
            })

        else:
            size = seg_end - seg_start + 1

            # Size check: blocks larger than gc_max cannot be gene conversions.
            # Classify by flanking context:
            #   het on both sides → DCO-internal (double crossover bracketing the region)
            #   LOH neighbor on either side → complex (ambiguous, examine by case)
            if size > gc_max:
                if fl in ("het", "tel") and fr in ("het", "tel"):
                    event_type = "DCO-internal"
                else:
                    event_type = "NCO-GC"  # placeholder; complex flag set below

                adj = _adjacent_to_terminal(fl, fr)

                rec = {
                    "type":                 event_type,
                    "chrom":                chrom["name"],
                    "start":                seg_start,
                    "end":                  seg_end,
                    "haplotype":            hap,
                    "flanking_left":        fl,
                    "flanking_right":       fr,
                    "adjacent_to_terminal": adj,
                }
                if event_type != "DCO-internal":
                    rec["complex"] = True
                events.append(rec)
                continue

            # Phase-aware CO/NCO detection for blocks within gc_max size.
            # Query homolog A haplotype just outside the LOH block on each side.
            # Use positions seg_start-1 and seg_end+1 (safe because block is
            # internal — neither boundary is 1 or chrom_length).
            hap_A_left  = _hap_A_at(seg_start - 1)
            hap_A_right = _hap_A_at(seg_end   + 1)

            # If hom A identity is the same on both sides: NCO (no phase switch).
            # If it differs: CO (phase switch — hom A and hom B swap across block).
            event_type = "NCO-GC" if hap_A_left == hap_A_right else "CO-GC"

            adj = _adjacent_to_terminal(fl, fr)

            events.append({
                "type":                 event_type,
                "chrom":                chrom["name"],
                "start":                seg_start,
                "end":                  seg_end,
                "haplotype":            hap,
                "flanking_left":        fl,
                "flanking_right":       fr,
                "adjacent_to_terminal": adj,
            })

    return events


def reclassify_terminal_clusters(events, gc_max=5000):
    """
    Apply heuristic reclassification to contiguous terminal LOH clusters.

    A terminal cluster is a maximal run of LOH blocks on one arm with no
    het gaps where at least one block touches a telomere.

    Rules
    -----
    Two-block cluster:
        Both blocks → CO-terminal. No further processing.

    Three-or-more block cluster:
        Interior blocks (flanked by LOH on both sides, not touching het/tel)
        are candidates for NCO-GC extraction if:
          - size <= gc_max, AND
          - both immediate LOH neighbors share the same haplotype

        If neighbors differ → block tagged "complex" (left as CO-terminal).

        After NCO-GC extraction, adjacent same-haplotype CO-terminal blocks
        are merged into a single record spanning their combined range.
        The merged record inherits the flanking values of the outermost
        original blocks (left flank of first, right flank of last).

        Blocks > gc_max or touching het remain CO-terminal individually.

    The adjacent_to_terminal flag is preserved on extracted NCO-GC records.
    All flanking values are from the original LOH map neighbors.

    Parameters
    ----------
    events : list of dict
        Output of classify_events().
    gc_max : int
        Maximum gene conversion tract size used in the simulation.

    Returns
    -------
    list of dict
        Reclassified events sorted by start position. New/modified events
        have a 'reclassified' key set to True. A 'complex' key is set to
        True on blocks that could not be cleanly resolved.
    """

    # -----------------------------------------------------------------------
    # Partition events into terminal clusters and non-cluster events
    # -----------------------------------------------------------------------
    # A terminal cluster is a maximal contiguous sequence of CO-terminal or
    # adjacent-to-terminal events with no het gap between them.
    # We identify them by grouping events whose flanking values connect
    # them to each other (either flanking_left or flanking_right is a
    # haplotype value, not "het" or "tel") and at least one touches a tel.

    # Work on a copy sorted by start position
    events = sorted(events, key=lambda e: e["start"])

    # Find indices of all events that are in a terminal cluster.
    # Strategy: collect runs of events where each is adjacent_to_terminal
    # or CO-terminal and they form a contiguous LOH run (no het between them).
    def _in_terminal_run(e):
        return e["type"] == "CO-terminal" or e.get("adjacent_to_terminal", False)

    # Group consecutive events that are part of terminal runs
    groups   = []   # list of (is_terminal_cluster, [event_indices])
    i        = 0
    n        = len(events)

    while i < n:
        if _in_terminal_run(events[i]):
            # Collect the full contiguous terminal run
            run = [i]
            j   = i + 1
            while j < n and _in_terminal_run(events[j]):
                run.append(j)
                j += 1
            # Only treat as a cluster if it actually contains a CO-terminal
            has_terminal = any(events[k]["type"] == "CO-terminal" for k in run)
            if has_terminal:
                groups.append(("cluster", run))
            else:
                for k in run:
                    groups.append(("solo", [k]))
            i = j
        else:
            groups.append(("solo", [i]))
            i += 1

    # -----------------------------------------------------------------------
    # Process each group
    # -----------------------------------------------------------------------
    result = []

    for kind, indices in groups:
        if kind == "solo":
            result.append(events[indices[0]])
            continue

        cluster = [events[k] for k in indices]
        m       = len(cluster)

        # Two-block rule: both become CO-terminal
        if m == 2:
            for e in cluster:
                rec = dict(e)
                rec["type"]         = "CO-terminal"
                rec["reclassified"] = True
                result.append(rec)
            continue

        # Three-or-more: extract NCO-GC islands then merge CO-terminal runs
        # ------------------------------------------------------------------
        # Pass 1: mark each block as NCO-GC, CO-terminal, or complex
        # ------------------------------------------------------------------
        labels = []   # "co-terminal" | "gc-nco" | "complex"
        for idx, e in enumerate(cluster):
            # Blocks touching het or tel on either side cannot be NCO-GC
            if e["flanking_left"] in ("het", "tel") or \
               e["flanking_right"] in ("het", "tel"):
                labels.append("co-terminal")
                continue

            # Interior block: check size and neighbor haplotypes
            size         = e["end"] - e["start"] + 1
            left_hap     = e["flanking_left"]   # haplotype of left LOH neighbor
            right_hap    = e["flanking_right"]  # haplotype of right LOH neighbor

            if size <= gc_max and left_hap == right_hap:
                labels.append("gc-nco")
            elif size <= gc_max and left_hap != right_hap:
                labels.append("complex")
            else:
                labels.append("co-terminal")

        # ------------------------------------------------------------------
        # Pass 2: merge adjacent co-terminal blocks of the same haplotype
        # ------------------------------------------------------------------
        # Build a list of (label, event) pairs after extraction
        tagged = list(zip(labels, cluster))

        # Emit NCO-GC and complex blocks as individual records;
        # collect runs of co-terminal blocks for merging
        co_run = []   # accumulator for current CO-terminal run

        def _flush_co_run(run):
            """Merge a run of CO-terminal blocks into one record."""
            if not run:
                return []
            if len(run) == 1:
                rec = dict(run[0])
                rec["type"]         = "CO-terminal"
                rec["reclassified"] = True
                return [rec]
            # All same haplotype (guaranteed by merge logic below)
            merged = dict(run[0])
            merged["type"]           = "CO-terminal"
            merged["end"]            = run[-1]["end"]
            merged["flanking_right"] = run[-1]["flanking_right"]
            merged["reclassified"]   = True
            return [merged]

        out_records = []
        for label, e in tagged:
            if label == "co-terminal":
                # Only merge with run if same haplotype as current run
                if co_run and co_run[-1]["haplotype"] != e["haplotype"]:
                    out_records.extend(_flush_co_run(co_run))
                    co_run = []
                co_run.append(e)
            else:
                # Flush any pending CO-terminal run first
                out_records.extend(_flush_co_run(co_run))
                co_run = []
                rec = dict(e)
                rec["type"]         = "NCO-GC" if label == "gc-nco" else e["type"]
                rec["reclassified"] = True
                if label == "complex":
                    rec["complex"]  = True
                    rec["type"]     = "CO-terminal"
                out_records.append(rec)

        out_records.extend(_flush_co_run(co_run))

        # Second pass: merge CO-terminal records of the same haplotype that
        # are separated only by NCO-GC islands (the islands are preserved).
        final_records = []
        for rec in out_records:
            if rec["type"] == "NCO-GC":
                # Check if the last CO-terminal and the next CO-terminal
                # (after this island) will be same haplotype — handled lazily:
                # just append and let the CO-terminal merge on the next CO hit
                final_records.append(rec)
            elif rec["type"] == "CO-terminal":
                # Find the most recent CO-terminal in final_records (skipping
                # any intervening NCO-GC islands)
                prev_co_idx = None
                for k in range(len(final_records) - 1, -1, -1):
                    if final_records[k]["type"] == "CO-terminal":
                        prev_co_idx = k
                        break
                    elif final_records[k]["type"] != "NCO-GC":
                        break  # hit something else — don't merge

                if (prev_co_idx is not None
                        and final_records[prev_co_idx]["haplotype"] == rec["haplotype"]):
                    # Merge: extend the previous CO-terminal
                    final_records[prev_co_idx]["end"]            = rec["end"]
                    final_records[prev_co_idx]["flanking_right"] = rec["flanking_right"]
                else:
                    final_records.append(rec)
            else:
                final_records.append(rec)

        result.extend(final_records)

    return sorted(result, key=lambda e: e["start"])


def _format_events(events):
    """
    Return a formatted string summary of classified events for CLI output.
    """
    if not events:
        return "  (no events detected)"
    lines = []
    for e in events:
        tags = []
        if e.get("adjacent_to_terminal"):
            tags.append("adjacent-to-terminal")
        if e.get("reclassified"):
            tags.append("reclassified")
        if e.get("complex"):
            tags.append("complex")
        tag_str = ("  [" + ", ".join(tags) + "]") if tags else ""
        lines.append(
            f"  {e['type']:<12}  "
            f"{e['start']:>8} – {e['end']:>8}  "
            f"[{e['haplotype']}]  "
            f"left={e['flanking_left']}  right={e['flanking_right']}"
            f"{tag_str}"
        )
    return "\n".join(lines)


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
    from datetime import datetime

    __version__ = "0.10"

    parser = argparse.ArgumentParser(description="Mitotic recombination simulator")
    parser.add_argument("genome_file", help="CSV with columns: chromosome, length, CEN")
    parser.add_argument("--version", action="version", version=f"simRec {__version__}")
    parser.add_argument("--n-gen", type=int, default=10, dest="n_gen", help="Number of generations")
    parser.add_argument("--p-rec", type=float, default=1e-5, dest="p_rec", help="Base recombination probability per bp")
    parser.add_argument("--gc-min", type=int, default=100, dest="gc_min", help="Min GC tract size (bp)")
    parser.add_argument("--gc-max", type=int, default=5000, dest="gc_max", help="Max GC tract size (bp)")
    parser.add_argument("--co-nco", type=float, default=0.5, dest="co_prob",
                        help="Probability of crossover per recombination event (default: 0.5, i.e. 1:1 CO:NCO)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--plot", action="store_true", help="Plot final haplotype and LOH state")
    parser.add_argument("--plot-out", type=str, default=None,
                        help="Save plot to this path (e.g. out.png) instead of displaying interactively")
    parser.add_argument("--observed", action="store_true",
                        help="Print post-hoc classified events (observed from the final LOH map) "
                             "to stdout as a tab-separated table with header")
    parser.add_argument("--observed-out", type=str, default=None, dest="observed_out",
                        metavar="FILE",
                        help="Write post-hoc classified events to FILE instead of stdout")
    parser.add_argument("--logged", action="store_true",
                        help="Print logged recombination events (tracked as the simulation "
                             "runs) to stdout as a tab-separated table with header")
    parser.add_argument("--logged-out", type=str, default=None, dest="logged_out",
                        metavar="FILE",
                        help="Write logged recombination events to FILE; also writes "
                             "haplotype.txt alongside. Columns: gen, chrom, type, start, end")
    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    genome = load_genome(args.genome_file, p_rec_default=args.p_rec)

    def _write_event_log(event_log, path):
        """Write the logged event log to a CSV file."""
        with open(path, "w") as fh:
            fh.write("gen,chrom,type,start,end,converts_to\n")
            for e in event_log:
                fh.write(f"{e['gen']},{e['chrom']},{e['type']},"
                         f"{e['start']},{e['end']},{e['converts_to']}\n")

    def _write_haplotype_map(haplotype_snapshots, log_path):
        """
        Write a text listing of haplotype segments for each generation in
        which events fired. Only generations with events are written.
        Written to haplotype.txt in the same directory as the log file.
        """
        import os
        hap_path = os.path.join(os.path.dirname(os.path.abspath(log_path)), "haplotype.txt")
        with open(hap_path, "w") as fh:
            for gen, cell in haplotype_snapshots:
                chrom_name = cell["A"]["name"]
                fh.write(f"Generation {gen}  Chr {chrom_name}\n")
                for homolog in ("A", "B"):
                    fh.write(f"  Homolog {homolog}:\n")
                    for seg_start, seg_end, hap in cell[homolog]["haplotype"]:
                        fh.write(f"    {seg_start:>8} – {seg_end:>8}  [{hap}]\n")
                fh.write("\n")

    def _format_logged_row(e):
        return "\t".join((
            str(e["gen"]), e["chrom"], e["type"],
            str(e["start"]), str(e["end"]), e["converts_to"],
        ))

    # Run simulation whenever any output mode is active, or in default verbose mode
    wants_observed = args.observed or args.observed_out
    wants_logged   = args.logged   or args.logged_out

    final_cell, event_log, haplotype_snapshots = run_simulation(
        genome, n_gen=args.n_gen, gc_min=args.gc_min,
        gc_max=args.gc_max, co_prob=args.co_prob)

    # --logged / --logged-out
    if args.logged:
        log_cols = ("gen", "chrom", "type", "start", "end", "converts_to")
        print("\t".join(log_cols))
        for e in event_log:
            print(_format_logged_row(e))

    if args.logged_out is not None:
        _write_event_log(event_log, args.logged_out)
        _write_haplotype_map(haplotype_snapshots, args.logged_out)

    # --observed / --observed-out
    if wants_observed:
        events    = classify_events(final_cell, gc_max=args.gc_max)
        reclass   = reclassify_terminal_clusters(events, gc_max=args.gc_max)
        timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        obs_cols  = ("time", "chrom", "event", "start", "end", "haplotype", "left", "right",
                     "adjacent_to_terminal", "reclassified", "complex")

        def _format_observed_row(e):
            return "\t".join((
                timestamp,
                e["chrom"],
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

        if args.observed:
            print("\t".join(obs_cols))
            for e in reclass:
                print(_format_observed_row(e))

        if args.observed_out is not None:
            with open(args.observed_out, "w") as fh:
                fh.write("\t".join(obs_cols) + "\n")
                for e in reclass:
                    fh.write(_format_observed_row(e) + "\n")

    # Default verbose mode — only when no TSV output flags are given
    if not wants_observed and not wants_logged:
        print(f"Loaded chromosome '{genome['A']['name']}', length {genome['A']['length']} bp")
        print(f"Running {args.n_gen} generation(s)...")

        print("\nFinal haplotype state:")
        for homolog in ("A", "B"):
            print(f"  Homolog {homolog}:")
            for seg in final_cell[homolog]["haplotype"]:
                print(f"    {seg[0]:>8} – {seg[1]:>8}  [{seg[2]}]")

        print("\nLOH state:")
        for seg in compute_loh(final_cell):
            print(f"    {seg[0]:>8} – {seg[1]:>8}  [{seg[2]}]")

        print("\nClassified events:")
        print(_format_events(classify_events(final_cell)))

        print("\nReclassified events:")
        print(_format_events(reclassify_terminal_clusters(
            classify_events(final_cell), gc_max=args.gc_max)))

    if args.plot or args.plot_out:
        plot_cell_and_loh(
            final_cell,
            title=f"After {args.n_gen} generation(s)",
            savepath=args.plot_out
        )
