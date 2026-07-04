"""Theoretical isotope distributions for peptides.

This is a faithful port of the validated ``individual_element_binomial_walk`` /
``descending_partial_products`` exploration (provided alongside the viewer spec)
onto the shared constants in ``chemistry.py``. The algorithm walks each
element's binomial isotope expansion as a priority heap, then takes a descending
partial-products merge across elements to enumerate the most-abundant
sub-formulas down to ``dividing_threshold`` of the base peak.

Public entry points:
  * ``peptide_isotope_distribution(seq)``    -> (neutral masses, abundances)
        per nominal-mass isotopologue group (the MS1 envelope).
  * ``peptide_isotope_mzs(seq, charge)``     -> (m/z, normalized abundances)
        ready to overlay against an experimental distribution in panel 3.

Only the public functions and the two ported workhorses are documented here;
the inner heap bookkeeping is unchanged from the original so behavior matches
the pipeline that produced the sqlite distributions.
"""

import bisect
import heapq
import itertools
import math
from collections import defaultdict

import numpy as np

try:
    from . import chemistry as chem
except ImportError:
    import chemistry as chem

# Aliases keep the ported algorithm identical to the original module-global form.
elementalprobabilities = chem.ELEMENTAL_PROBABILITIES
elementalmasses = chem.ELEMENTAL_MASSES
isotopesbyelement = chem.ISOTOPES_BY_ELEMENT
monoisotopickeys = chem.MONOISOTOPIC_KEYS
nonmonoisotopicgroups = chem.NONMONOISOTOPIC_GROUPS
nonmonoisotopicelements = chem.NONMONOISOTOPIC_ELEMENTS
elementvectors = chem.ELEMENT_VECTORS
vectorpositions = chem.VECTOR_POSITIONS
elementpositions = chem.ELEMENT_POSITIONS
proton = chem.PROTON

DEFAULT_DIVIDING_THRESHOLD = 0.01
DEFAULT_SUBISO_DEPTH = 0.8


def individual_element_binomial_walk(dividingthreshold, e, acount):
    """Enumerate the isotopologues of ``acount`` atoms of element ``e``.

    Returns a heapified list of ``[ratio, prob, mass, element, vector]`` entries
    covering the abundant isotope combinations down to the threshold.
    """
    elementlist = []
    mainheap = []
    vectorsets = defaultdict(set)
    mk = monoisotopickeys[e]
    nvector = elementvectors[e].copy()
    nvector[vectorpositions[e][mk]] += acount
    if len(isotopesbyelement[e]) > 2:
        baseprob = elementalprobabilities[mk] ** acount
        preheap = []
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        while greater:
            greater = False
            for iso in nonmonoisotopicgroups[e]:
                newelementvector = nvector.copy()
                newelementvector[vectorpositions[e][mk]] -= 1
                if newelementvector[vectorpositions[e][mk]] > -1:
                    newelementvector[vectorpositions[e][iso]] += 1
                    vectorsets[e].add(tuple(newelementvector))
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(newelementvector):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso] ** c
                        if loopiso in nonmonoisotopicelements:
                            newelementprob *= math.comb(acount - pn, c)
                            pn += c
                    preheap.append([newelementprob, newelementmass, e, newelementvector.copy()])
                    if newelementprob > lastprob:
                        lastprob = newelementprob
                        greater = True
        preheap = sorted(preheap)
        maxiso = preheap[-1]
        maxprob, m, e, nv = maxiso
        elementlist.append([-1, maxprob, m, e, nv])
        maxprob *= -1
        preheap = preheap[:-1]
        for h in preheap:
            r = h[0] / maxprob
            h.insert(0, r)
            heapq.heappush(mainheap, h)
        for iso in nonmonoisotopicgroups[e]:
            v = nv.copy()
            v[vectorpositions[e][mk]] -= 1
            if v[vectorpositions[e][mk]] > -1:
                v[vectorpositions[e][iso]] += 1
                tuplevec = tuple(v)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(v):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso] ** c
                        if loopiso in nonmonoisotopicelements:
                            newelementprob *= math.comb(acount - pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    else:
        preheap = []
        baseprob = elementalprobabilities[mk] ** acount
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        iso = nonmonoisotopicgroups[e][0]
        while greater:
            greater = False
            nvector[vectorpositions[e][mk]] -= 1
            if nvector[vectorpositions[e][mk]] > -1:
                nvector[vectorpositions[e][iso]] += 1
                vectorsets[e].add(tuple(nvector))
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(nvector):
                    loopiso = elementpositions[e][n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= elementalprobabilities[loopiso] ** c
                    if loopiso in nonmonoisotopicelements:
                        newelementprob *= math.comb(acount - pn, c)
                        pn += c
                preheap.append([newelementprob, newelementmass, e, nvector.copy()])
                if newelementprob > lastprob:
                    lastprob = newelementprob
                    greater = True
        preheap = sorted(preheap)
        maxiso = preheap[-1]
        maxprob, m, e, nv = maxiso
        elementlist.append([-1, maxprob, m, e, nv])
        maxprob *= -1
        preheap = preheap[:-1]
        for h in preheap:
            r = h[0] / maxprob
            h.insert(0, r)
            heapq.heappush(mainheap, h)
        v = nv.copy()
        v[vectorpositions[e][mk]] -= 1
        if v[vectorpositions[e][mk]] > -1:
            v[vectorpositions[e][iso]] += 1
            tuplevec = tuple(v)
            if tuplevec not in vectorsets[e]:
                vectorsets[e].add(tuplevec)
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(v):
                    loopiso = elementpositions[e][n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= elementalprobabilities[loopiso] ** c
                    if loopiso in nonmonoisotopicelements:
                        newelementprob *= math.comb(acount - pn, c)
                        pn += c
                heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])

    cutoff = -maxprob * dividingthreshold

    r, p, m, e, v = heapq.heappop(mainheap)
    elementlist.append([r, p, m, e, v])
    if len(isotopesbyelement[e]) > 2:
        while p > cutoff:
            for iso in nonmonoisotopicgroups[e]:
                newelementvector = v.copy()
                newelementvector[vectorpositions[e][mk]] -= 1
                if newelementvector[vectorpositions[e][mk]] > 0:
                    newelementvector[vectorpositions[e][iso]] += 1
                    tuplevec = tuple(newelementvector)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(newelementvector):
                            loopiso = elementpositions[e][n]
                            newelementmass += elementalmasses[loopiso] * c
                            newelementprob *= elementalprobabilities[loopiso] ** c
                            if loopiso in nonmonoisotopicelements:
                                newelementprob *= math.comb(acount - pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, newelementvector.copy()])
            r, p, m, e, v = heapq.heappop(mainheap)
            elementlist.append([r, p, m, e, v])
    else:
        iso = nonmonoisotopicgroups[e][0]
        while p > cutoff:
            nvector = v.copy()
            nvector[vectorpositions[e][mk]] -= 1
            if nvector[vectorpositions[e][mk]] > 0:
                nvector[vectorpositions[e][iso]] += 1
                tuplevec = tuple(nvector)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(nvector):
                        loopiso = elementpositions[e][n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= elementalprobabilities[loopiso] ** c
                        if loopiso in nonmonoisotopicelements:
                            newelementprob *= math.comb(acount - pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, nvector.copy()])
            r, p, m, e, v = heapq.heappop(mainheap)
            elementlist.append([r, p, m, e, v])
    heapq.heapify(elementlist)
    return elementlist


def descending_partial_products(dividingthreshold, elementalorganizer):
    """Merge per-element isotopologue heaps into whole-molecule sub-formulas."""
    for k in elementalorganizer:
        heapq.heapify(elementalorganizer[k])

    mainpool = defaultdict(list)
    for k in elementalorganizer:
        mainpool[k].append(heapq.heappop(elementalorganizer[k]))

    formula = ""
    maxprob = 1
    mainmass = 0
    finalabundances = {}
    for b in sorted(mainpool):
        for r, p, m, e, v in mainpool[b]:
            for n, c in enumerate(v):
                if c > 0:
                    formula += f"{elementpositions[e][n]}({c})"
            maxprob *= p
            mainmass += m

    finalabundances[formula] = [mainmass, maxprob]

    cutoff = maxprob * dividingthreshold
    mainheap = list(itertools.chain(*elementalorganizer.values()))
    heapq.heapify(mainheap)

    multinomialpath = []
    probabilityranking = []
    while mainheap:
        r, p, m, e, v = heapq.heappop(mainheap)
        baseiter = {k: v for k, v in mainpool.items() if k != e}
        baseiter[e] = [[r, p, m, e, v]]

        formula = ""
        prob = 1
        mass = 0
        for b in sorted(baseiter):
            for sr, sp, sm, se, sv in baseiter[b]:
                for n, c in enumerate(sv):
                    if c > 0:
                        formula += f"{elementpositions[se][n]}({c})"
                prob *= sp
                mass += sm

        finalabundances[formula] = [mass, prob]
        if prob < cutoff:
            break

        ind = bisect.bisect(probabilityranking, r)
        probabilityranking.insert(ind, r)
        multinomialpath.insert(ind, [r, p, m, e, v])

        checkedcombos = set()
        for path in multinomialpath.copy():
            multielement = False
            match path[1]:
                case list():
                    multielement = True
                    sepool = set()
                    sepool.add(e)
                    seformulas = []
                    multipath = []
                    nsr = 1
                    for sr, sp, sm, se, sv in path[1:]:
                        if se not in sepool:
                            nsr *= sr
                            sepool.add(se)
                            sef = ""
                            for n, c in enumerate(sv):
                                if c > 0:
                                    sef += f"{elementpositions[se][n]}({c})"
                            seformulas.append(sef)
                            multipath.append([sr, sp, sm, se, sv])
                    checkformula = "".join(sorted(seformulas))
                    if checkformula in checkedcombos:
                        continue
                    else:
                        checkedcombos.add(checkformula)
                    if len(multipath) == 0:
                        continue
                case _:
                    sr, sp, sm, se, sv = path
                    sef = "".join(f"{se}{str(n)}{(val)}" for n, val in enumerate(sv))
                    if sef in checkedcombos:
                        continue
                    else:
                        checkedcombos.add(sef)
                    if se == e:
                        continue
                    nsr = sr
            newratio = nsr * r
            if newratio > 0:
                newratio *= -1
            if -newratio >= dividingthreshold:
                if multielement:
                    seformula = ""
                    newprob = 1
                    newmass = 0
                    newiter = {k: v for k, v in baseiter.items() if k not in sepool}
                    newiter[e] = [[r, p, m, e, v]]
                    for ir, ip, im, ie, iv in multipath:
                        newiter[ie] = [[ir, ip, im, ie, iv]]
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    seformula += f"{elementpositions[ie][n]}({c})"
                            newprob *= ip
                            newmass += im
                else:
                    newiter = {k: v for k, v in baseiter.items() if k != se}
                    newiter[se] = [[sr, sp, sm, se, sv]]
                    seformula = ""
                    newprob = 1
                    newmass = 0
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    seformula += f"{elementpositions[ie][n]}({c})"
                            newprob *= ip
                            newmass += im
                if newprob >= cutoff:
                    finalabundances[seformula] = [newmass, newprob]
                    if multielement:
                        ind = bisect.bisect(probabilityranking, newratio)
                        probabilityranking.insert(ind, newratio)
                        multinomialpath.insert(ind, [newratio, *multipath])
                    else:
                        ind = bisect.bisect(probabilityranking, newratio)
                        probabilityranking.insert(ind, newratio)
                        multinomialpath.insert(ind, [newratio, [sr, sp, sm, se, sv], [r, p, m, e, v]])
            else:
                break

    subformulas, massesandabundances = list(zip(*finalabundances.items()))
    subformulas = np.array(subformulas, dtype="S")
    massesandabundances = np.array(massesandabundances).transpose()
    subformulas = subformulas[massesandabundances[0].argsort()].tolist()
    massesandabundances = massesandabundances[:, massesandabundances[0].argsort()]
    return subformulas, massesandabundances


def distribution_generation(dividingthreshold, atomiccomposition):
    elementalorganizer = {}
    for e, acount in atomiccomposition.items():
        elementalorganizer[e] = individual_element_binomial_walk(dividingthreshold, e, acount).copy()
    return descending_partial_products(dividingthreshold, elementalorganizer)


def peptide_isotope_distribution(sequence, dividing_threshold=DEFAULT_DIVIDING_THRESHOLD,
                                 subiso_depth=DEFAULT_SUBISO_DEPTH):
    """Neutral isotope envelope of a peptide.

    Returns ``(masses, abundances)`` where each entry is one nominal-mass
    isotopologue group: ``masses`` are intensity-weighted mean neutral masses
    and ``abundances`` are the summed group abundances (not normalized).
    """
    atomiccomposition = chem.peptide_atomic_composition(sequence)
    subformulas, massesandabundances = distribution_generation(dividing_threshold, atomiccomposition)
    subformulas = [i.decode() for i in subformulas]

    massgroups = defaultdict(list)
    intensitygroups = defaultdict(list)
    masses, abundances = massesandabundances
    for n, s in enumerate(subformulas):
        massnumber = 0
        for ss in s.split(")")[:-1]:
            i1, i2 = map(int, ss[1:].split("("))
            massnumber += i1 * i2
        massgroups[massnumber].append(masses[n])
        intensitygroups[massnumber].append(abundances[n])

    meansofmasses = []
    sumsofabundances = []
    for mn, m in massgroups.items():
        a = intensitygroups[mn]
        totalabundance = sum(a)
        weightedmass = sum(sm * sa for sm, sa in zip(m, a))
        meansofmasses.append(weightedmass / totalabundance)
        sumsofabundances.append(totalabundance)

    order = np.argsort(meansofmasses)
    masses = np.array(meansofmasses)[order]
    abundances = np.array(sumsofabundances)[order]
    return masses, abundances


def peptide_isotope_mzs(sequence, charge, dividing_threshold=DEFAULT_DIVIDING_THRESHOLD):
    """Charged m/z and max-normalized abundances for overlay plots.

    The abundances are scaled so the most abundant isotopologue is 1.0, which is
    what panel 3 needs to twin-plot against an experimental peak at matched height.
    """
    charge = max(1, int(charge))
    masses, abundances = peptide_isotope_distribution(sequence, dividing_threshold)
    mzs = (masses + proton * charge) / charge
    peak = abundances.max() if abundances.size else 1.0
    return mzs, abundances / peak


def peptide_isotope_raw(sequence, dividing_threshold=DEFAULT_DIVIDING_THRESHOLD):
    """RAW isotopologue distribution: one entry per sub-formula, with NO joining
    of isotopomers that share a nominal M+N mass (cf. libraryadditions.py's raw
    ``massesandabundances`` before ``subisotopomer_handler`` groups them).

    Returns ``(masses, abundances)`` sorted by mass.
    """
    atomiccomposition = chem.peptide_atomic_composition(sequence)
    _subformulas, massesandabundances = distribution_generation(
        dividing_threshold, atomiccomposition)
    masses = np.asarray(massesandabundances[0], dtype=float)
    abundances = np.asarray(massesandabundances[1], dtype=float)
    order = masses.argsort()
    return masses[order], abundances[order]


def peptide_isotope_bars(sequence, charge, mode="raw",
                         dividing_threshold=DEFAULT_DIVIDING_THRESHOLD):
    """Theoretical MS1 bars (m/z, abundance) for an overlay.

    ``mode='raw'`` keeps every isotopomer separate; ``mode='summed'`` merges all
    isotopomers at the same M+N into one signal (the ``sumabundancedist``). m/z
    is at the given precursor ``charge``.
    """
    charge = max(1, int(charge or 1))
    if mode == "summed":
        masses, abundances = peptide_isotope_distribution(sequence, dividing_threshold)
    else:
        masses, abundances = peptide_isotope_raw(sequence, dividing_threshold)
    mzs = (masses + proton * charge) / charge
    return mzs, abundances
