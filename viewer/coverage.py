"""Peptide sequence coverage by matched MS2 fragment ions.

A faithful port of the coverage concept in
``examples/sequencecoverageconcept.py`` (the ``coverage_print`` /
divider-string logic), adapted for the viewer's Table 2 so each candidate
peptide can be shown with the residue coverage its fragment ions provide.

The idea: every matched b/y ion "cuts" the sequence at a backbone position.
N-terminal ions (a/b/c) cover the residues left of the cut; C-terminal ions
(x/y/z) cover the residues to the right. Collecting all the cut points splits
the sequence into segments (the divider string, ``A|BC|DEF``), and each segment
is "covered" by however many ions span it. A peptide whose ions isolate small,
heavily-covered segments is more confidently distinguishable from a decoy than
one whose ions only chop it into a couple of big pieces.

Public entry points:
  * ``fragment_ions(seq, charges, ions)``  -> {ion_label: m/z} theoretical b/y.
  * ``match_coverage(seq, mz, ..., ppm)``  -> the list of ion labels (``b3`` …)
        whose theoretical m/z matched a peak in the MS2 spectrum.
  * ``coverage_string(seq, ion_labels)``   -> the ``A|BC|DEF`` divider string.
  * ``coverage_summary(seq, mz, ...)``     -> a one-line cell string for Table 2.
"""

from collections import Counter, defaultdict

import numpy as np

try:
    from . import chemistry as chem
except ImportError:
    import chemistry as chem

proton = chem.PROTON

DEFAULT_IONS = "by"
DEFAULT_PPM = 20.0
DEFAULT_CHARGES = (1, 2)


def fragment_ions(seq, charges=DEFAULT_CHARGES, ions=DEFAULT_IONS):
    """Theoretical fragment-ion m/z for a bare peptide ``seq``.

    Returns ``{ion_label: [m/z per charge]}`` where ``ion_label`` is the ion
    type + backbone position (``b3``, ``y4``). N-terminal ions accumulate
    residue compositions from the start, C-terminal ions from the end; the
    per-ion-type composition deltas come from ``chemistry`` (same values the
    reference used). Mirrors ``fragmentation_compositions`` in the reference.
    """
    seq = seq or ""
    n = len(seq)
    out = {}
    if n < 2:
        return out

    n_ions = [i for i in ions if i in chem.N_FRAGMENT_COMPOSITIONS]
    c_ions = [i for i in ions if i in chem.C_FRAGMENT_COMPOSITIONS]

    # N-terminal: cumulative residues seq[:k] for k = 1 .. n-1.
    running = Counter()
    for k in range(1, n):
        running += Counter(chem.AMINO_ACID_COMPOSITION[seq[k - 1]])
        for ion in n_ions:
            comp = running + Counter(chem.N_FRAGMENT_COMPOSITIONS[ion])
            mass = chem.monoisotopic_mass(comp)
            out[f"{ion}{k}"] = [(mass + proton * z) / z for z in charges]

    # C-terminal: cumulative residues from the end, seq[-k:] for k = 1 .. n-1.
    running = Counter()
    for k in range(1, n):
        running += Counter(chem.AMINO_ACID_COMPOSITION[seq[n - k]])
        for ion in c_ions:
            comp = running + Counter(chem.C_FRAGMENT_COMPOSITIONS[ion])
            mass = chem.monoisotopic_mass(comp)
            out[f"{ion}{k}"] = [(mass + proton * z) / z for z in charges]

    return out


def match_coverage(seq, peak_mz, charges=DEFAULT_CHARGES, ions=DEFAULT_IONS,
                   ppm=DEFAULT_PPM):
    """Ion labels whose theoretical m/z matched a peak in ``peak_mz``.

    ``peak_mz`` is the array of MS2 spectrum peak m/z. An ion counts as matched
    if any of its charge states falls within ``ppm`` of a peak. Returns a
    de-duplicated, sequence-ordered list of ``(type, position)`` labels
    (``b3`` …) — the input to :func:`coverage_string`.
    """
    peaks = np.asarray(peak_mz, dtype=float)
    if peaks.size == 0:
        return []
    peaks = np.sort(peaks)

    matched = []
    for label, mzs in fragment_ions(seq, charges, ions).items():
        for mz in mzs:
            tol = mz / 1e6 * ppm
            i = int(np.searchsorted(peaks, mz))
            near = False
            for j in (i - 1, i):
                if 0 <= j < peaks.size and abs(peaks[j] - mz) <= tol:
                    near = True
                    break
            if near:
                matched.append(label)
                break
    return matched


def coverage_string(seq, ion_labels):
    """Render the ``A|BC|DEF`` divider string for matched ion labels.

    Ported from ``coverage_print``: each N-terminal ion cuts after ``ioncount``
    residues, each C-terminal ion cuts ``len-ioncount`` from the start; the set
    of cuts splits the sequence into segments joined by ``|``.
    """
    olen = len(seq)
    dividers = []
    for ion in ion_labels:
        iontype = ion[0]
        ioncount = int(ion[1:])
        if iontype in "abc":      # N-terminal: residues [0:ioncount]
            dividers.append(ioncount)
        elif iontype in "xyz":    # C-terminal: residues [olen-ioncount:]
            dividers.append(olen - ioncount)
    dividers = sorted(set(d for d in dividers if 0 < d < olen))
    if not dividers:
        return seq

    parts = []
    prev = 0
    for d in dividers:
        parts.append(seq[prev:d])
        prev = d
    parts.append(seq[prev:])
    return "|".join(p for p in parts if p)


def coverage_segments(seq, ion_labels):
    """Per-segment coverage counts: ``[(segment, n_ions_spanning), ...]``.

    A segment is covered by an N-terminal ion whose cut lies to its right and a
    C-terminal ion whose cut lies at or before its left edge (the reference's
    ``ntermcovers``/``ctermcovers`` test). Useful for shading/tooltips beyond
    the flat divider string.
    """
    olen = len(seq)
    nterm = []  # cut positions of N-terminal ions
    cterm = []  # cut positions of C-terminal ions
    dividers = []
    for ion in ion_labels:
        iontype = ion[0]
        ioncount = int(ion[1:])
        if iontype in "abc":
            dividers.append(ioncount)
            nterm.append(ioncount)
        elif iontype in "xyz":
            dividers.append(olen - ioncount)
            cterm.append(olen - ioncount)
    dividers = sorted(set(d for d in dividers if 0 < d < olen))

    segments = []
    bounds = [0] + dividers + [olen]
    for a, b in zip(bounds[:-1], bounds[1:]):
        if a == b:
            continue
        covers = sum(1 for i in nterm if i > a) + sum(1 for i in cterm if i <= a)
        segments.append((seq[a:b], covers))
    return segments


def coverage_metrics(seq, ion_labels):
    """The coverage SCORE and its components, ported verbatim from the metric
    block in ``sequencecoverageconcept.py`` (the ``finalmetrics`` computation).

    The fragmentation isolates the sequence into segments (``partialseqs``);
    each segment that is spanned by at least one ion is "covered" ``count``
    times. The reference scores a peptide by how finely and redundantly its
    ions isolate the backbone:

      * ``coverageweight``       = 1 / (max N-term cut + max C-term cut)
      * ``isolationlengthweight``= ∏ over covered segments of 1/len(seg)/len(seq)
      * ``dividerweight``        = 1 / (number of covered segments)
      * ``score``                = dividerweight·isolationlengthweight·coverageweight

    Returns a dict with ``score`` and the three component weights plus the
    segment count and a ``matchcounts`` tally (``secondfinalmetrics``). Returns
    all-zero metrics when nothing is covered.
    """
    slen = len(seq)
    zero = {"score": 0.0, "dividerweight": 0.0, "isolationlengthweight": 0.0,
            "coverageweight": 0.0, "nsegments": 0, "matchcounts": 0}
    if slen < 1 or not ion_labels:
        return zero

    maxncoverage = 0   # max C-term cut (reference's confusingly-named var)
    maxccoverage = 0   # max N-term cut
    dividers = []
    ntermcoverage = []
    ctermcoverage = []
    for ion in ion_labels:
        iontype = ion[0]
        ioncount = int(ion[1:])
        if iontype in "abc":
            dividers.append(ioncount)
            ntermcoverage.append(ioncount)
            if ioncount > maxccoverage:
                maxccoverage = ioncount
        elif iontype in "xyz":
            dividers.append(slen - ioncount)
            ctermcoverage.append(slen - ioncount)
            if ioncount > maxncoverage:
                maxncoverage = ioncount
    if (maxncoverage + maxccoverage) == 0:
        return zero
    dividers = sorted(set(dividers))
    coverageweight = 1.0 / (maxncoverage + maxccoverage)

    # Isolated segments (start..start+d), keyed by start index so identical
    # partial sequences at different positions stay distinct.
    ind = 0
    ddiff = np.diff(dividers, prepend=0).tolist()
    partialseqs = defaultdict(int)
    for d in ddiff:
        pseq = seq[ind:ind + d]
        ntermcovers = [i for i in ntermcoverage if i > ind]
        ctermcovers = [i for i in ctermcoverage if i <= ind]
        covers = len(ntermcovers) + len(ctermcovers)
        if covers > 0:
            partialseqs[f"{ind}-{pseq}"] += covers
        ind += d
    pseq = seq[ind:]
    ntermcovers = [i for i in ntermcoverage if i > ind]
    ctermcovers = [i for i in ctermcoverage if i <= ind]
    covers = len(ntermcovers) + len(ctermcovers)
    if covers > 0:
        partialseqs[f"{ind}-{pseq}"] += covers

    if not partialseqs:
        return zero

    matchcounts = len(set(ion_labels))
    isolationlengthweight = 1.0
    for indseq, count in partialseqs.items():
        _ind, pseq = indseq.split("-", 1)
        isolationlengthweight *= 1.0 / len(pseq) / len(seq)
        matchcounts += len(pseq) * count
    dividerweight = 1.0 / len(partialseqs)
    score = dividerweight * isolationlengthweight * coverageweight
    return {"score": score, "dividerweight": dividerweight,
            "isolationlengthweight": isolationlengthweight,
            "coverageweight": coverageweight,
            "nsegments": len(partialseqs), "matchcounts": matchcounts}


def coverage_summary(seq, peak_mz, charges=DEFAULT_CHARGES, ions=DEFAULT_IONS,
                     ppm=DEFAULT_PPM):
    """One-line Table-2 cell for a candidate peptide.

    Combines the matched-ion count, the residue coverage fraction, and the
    divider string, e.g. ``5 ions  6/8 aa  PEP|TIDE``. Returns ``"no match"``
    when nothing matched so empty candidates read clearly.
    """
    seq = seq or ""
    if len(seq) < 2:
        return ""
    matched = match_coverage(seq, peak_mz, charges, ions, ppm)
    if not matched:
        return "no match"
    covered = _covered_residues(seq, matched)
    return (f"{len(matched)} ions  {covered}/{len(seq)} aa  "
            f"{coverage_string(seq, matched)}")


def coverage_report(seq, peak_mz, charges=DEFAULT_CHARGES, ions=DEFAULT_IONS,
                    ppm=DEFAULT_PPM):
    """Everything Table 2 needs for one candidate peptide, in a single pass.

    Returns a dict: ``matched`` (ion labels), ``divider`` (the ``PEP|T|I|DER``
    string), ``covered`` (residues touched), ``score`` and the metric
    components from :func:`coverage_metrics`. ``score`` is 0 with an empty
    divider when nothing matched.
    """
    seq = seq or ""
    matched = match_coverage(seq, peak_mz, charges, ions, ppm) if len(seq) >= 2 else []
    metrics = coverage_metrics(seq, matched)
    return {
        "matched": matched,
        "divider": coverage_string(seq, matched) if matched else "",
        "covered": _covered_residues(seq, matched),
        **metrics,
    }


def _covered_residues(seq, ion_labels):
    """Number of residues touched by at least one matched ion."""
    olen = len(seq)
    covered = set()
    for ion in ion_labels:
        iontype = ion[0]
        ioncount = int(ion[1:])
        if iontype in "abc":
            covered.update(range(0, min(ioncount, olen)))
        elif iontype in "xyz":
            covered.update(range(max(0, olen - ioncount), olen))
    return len(covered)
