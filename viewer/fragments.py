"""Theoretical MS2 fragment-ion isotope masses for a peptide.

A faithful port of the fragment-isotope machinery from
``examples/peptidefragmentscoring.py`` /
``examples/sequencecoverageconcept.py``: ``fragment_element_binomial_walk`` +
``fragment_descending_partial_products`` enumerate a fragment ion's isotopomers
*constrained by the parent precursor's isotope composition*, so each fragment
mass can be attributed to the MS1 (precursor) isotopologue it descended from.

High-level entry point:
  * ``peptide_fragment_ions(seq, charge, iso_low, iso_high, ...)`` ->
        ``[{ion, isotope_index, neutral_mass, mz}]`` for the b/y ions of every
        precursor isotope whose m/z falls inside the MS2 isolation window
        ``[iso_low, iso_high]``. ``isotope_index`` is the precursor isotope the
        fragment came from (0 = monoisotopic, 1 = M+1, ...).

The heap bookkeeping in the two ported walkers is unchanged from the original so
behaviour matches the validated reference; only constant access is rewired onto
``chemistry.py``.
"""

import bisect
import heapq
import itertools
import math
from collections import defaultdict

import numpy as np

try:
    from . import chemistry as chem
    from .isotopes import distribution_generation
except ImportError:
    import chemistry as chem
    from isotopes import distribution_generation

elementalmasses = chem.ELEMENTAL_MASSES
proton = chem.PROTON

DEFAULT_DIVIDING_THRESHOLD = 0.1
DEFAULT_SUBISO_DEPTH = 0.5


# ---- ported fragment isotope walkers (verbatim logic) ----------------------

def fragment_element_binomial_walk(dividingthreshold, e, acount, fragprobabilities):
    nvector = []
    fragmentvectorpositions = {}
    fragmentelementpositions = {}
    maxinitial = 0
    mk = None
    for n, (iso, prob) in enumerate(fragprobabilities.items()):
        nvector.append(0)
        fragmentvectorpositions[iso] = n
        fragmentelementpositions[n] = iso
        if prob > maxinitial:
            maxinitial = prob
            mk = iso
    lesserfragmentisotopes = [i for i in fragprobabilities if i != mk]
    elementlist = []
    mainheap = []
    vectorsets = defaultdict(set)
    nvector[fragmentvectorpositions[mk]] += acount
    flen = len(fragprobabilities)
    if flen > 2:
        baseprob = fragprobabilities[mk] ** acount
        preheap = []
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        while greater:
            greater = False
            for iso in lesserfragmentisotopes:
                newelementvector = nvector.copy()
                newelementvector[fragmentvectorpositions[mk]] -= 1
                if newelementvector[fragmentvectorpositions[mk]] > -1:
                    newelementvector[fragmentvectorpositions[iso]] += 1
                    vectorsets[e].add(tuple(newelementvector))
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(newelementvector):
                        loopiso = fragmentelementpositions[n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso] ** c
                        if n > 0:
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
        for iso in lesserfragmentisotopes:
            v = nv.copy()
            v[fragmentvectorpositions[mk]] -= 1
            if v[fragmentvectorpositions[mk]] > -1:
                v[fragmentvectorpositions[iso]] += 1
                tuplevec = tuple(v)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(v):
                        loopiso = fragmentelementpositions[n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso] ** c
                        if n > 0:
                            newelementprob *= math.comb(acount - pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])
    else:
        preheap = []
        baseprob = fragprobabilities[mk] ** acount
        preheap.append([baseprob, acount * elementalmasses[mk], e, nvector.copy()])
        greater = True
        lastprob = baseprob
        iso = lesserfragmentisotopes[0]
        while greater:
            greater = False
            nvector[fragmentvectorpositions[mk]] -= 1
            if nvector[fragmentvectorpositions[mk]] > -1:
                nvector[fragmentvectorpositions[iso]] += 1
                vectorsets[e].add(tuple(nvector))
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(nvector):
                    loopiso = fragmentelementpositions[n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= fragprobabilities[loopiso] ** c
                    if n > 0:
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
        v[fragmentvectorpositions[mk]] -= 1
        if v[fragmentvectorpositions[mk]] > -1:
            v[fragmentvectorpositions[iso]] += 1
            tuplevec = tuple(v)
            if tuplevec not in vectorsets[e]:
                vectorsets[e].add(tuplevec)
                pn = 0
                newelementmass = 0
                newelementprob = 1
                for n, c in enumerate(v):
                    loopiso = fragmentelementpositions[n]
                    newelementmass += elementalmasses[loopiso] * c
                    newelementprob *= fragprobabilities[loopiso] ** c
                    if n > 0:
                        newelementprob *= math.comb(acount - pn, c)
                        pn += c
                heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, v.copy()])

    cutoff = -maxprob * dividingthreshold

    r, p, m, e, v = heapq.heappop(mainheap)
    elementlist.append([r, p, m, e, v])
    if flen > 2:
        while p > cutoff:
            for iso in lesserfragmentisotopes:
                newelementvector = v.copy()
                newelementvector[fragmentvectorpositions[mk]] -= 1
                if newelementvector[fragmentvectorpositions[mk]] > 0:
                    newelementvector[fragmentvectorpositions[iso]] += 1
                    tuplevec = tuple(newelementvector)
                    if tuplevec not in vectorsets[e]:
                        vectorsets[e].add(tuplevec)
                        pn = 0
                        newelementmass = 0
                        newelementprob = 1
                        for n, c in enumerate(newelementvector):
                            loopiso = fragmentelementpositions[n]
                            newelementmass += elementalmasses[loopiso] * c
                            newelementprob *= fragprobabilities[loopiso] ** c
                            if n > 0:
                                newelementprob *= math.comb(acount - pn, c)
                                pn += c
                        heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, newelementvector.copy()])
            r, p, m, e, v = heapq.heappop(mainheap)
            elementlist.append([r, p, m, e, v])
            try:
                r, p, m, e, v = heapq.heappop(mainheap)
                elementlist.append([r, p, m, e, v])
            except IndexError:
                break
    else:
        iso = lesserfragmentisotopes[0]
        while p > cutoff:
            nvector = v.copy()
            nvector[fragmentvectorpositions[mk]] -= 1
            if nvector[fragmentvectorpositions[mk]] > 0:
                nvector[fragmentvectorpositions[iso]] += 1
                tuplevec = tuple(nvector)
                if tuplevec not in vectorsets[e]:
                    vectorsets[e].add(tuplevec)
                    pn = 0
                    newelementmass = 0
                    newelementprob = 1
                    for n, c in enumerate(nvector):
                        loopiso = fragmentelementpositions[n]
                        newelementmass += elementalmasses[loopiso] * c
                        newelementprob *= fragprobabilities[loopiso] ** c
                        if n > 0:
                            newelementprob *= math.comb(acount - pn, c)
                            pn += c
                    heapq.heappush(mainheap, [newelementprob / maxprob, newelementprob, newelementmass, e, nvector.copy()])
            try:
                r, p, m, e, v = heapq.heappop(mainheap)
                elementlist.append([r, p, m, e, v])
            except IndexError:
                break
    heapq.heapify(elementlist)
    return elementlist, fragmentelementpositions


def fragment_descending_partial_products(dividingthreshold, elementalorganizer, fragmentpositions):
    mainpool = defaultdict(list)
    for k in elementalorganizer:
        mainpool[k].append(heapq.heappop(elementalorganizer[k]))

    subformulas = []
    sumabundances = []
    massnumberindices = {}

    formula = ''
    maxprob = 1
    mainmass = 0
    massnumber = 0
    for b in sorted(mainpool):
        for r, p, m, e, v in mainpool[b]:
            for n, c in enumerate(v):
                if c > 0:
                    iso = fragmentpositions[e][n]
                    massnumber += int(iso[1:]) * c
                    formula += f'{iso}({c})'
            maxprob *= p
            mainmass += m

    massnumberindices[massnumber] = 0
    subformulas.append(formula)
    sumabundances.append([mainmass * maxprob, maxprob])

    cutoff = maxprob * dividingthreshold
    mainheap = list(itertools.chain(*elementalorganizer.values()))
    heapq.heapify(mainheap)

    vectorpool = set()
    multinomialpath = []
    probabilityranking = []
    while mainheap:
        r, p, m, e, v = heapq.heappop(mainheap)
        baseiter = {k: v for k, v in mainpool.items() if k != e}
        baseiter[e] = [(r, p, m, e, v)]

        formula = ''
        prob = 1
        mass = 0
        massnumber = 0
        for b in sorted(baseiter):
            for sr, sp, sm, se, sv in baseiter[b]:
                for n, c in enumerate(sv):
                    if c > 0:
                        iso = fragmentpositions[se][n]
                        massnumber += int(iso[1:]) * c
                        formula += f'{iso}({c})'
                prob *= sp
                mass += sm

        try:
            index = massnumberindices[massnumber]
            subformulas[index] += '-' + formula
            sumabundances[index][0] += mass * prob
            sumabundances[index][1] += prob
        except KeyError:
            index = len(massnumberindices)
            massnumberindices[massnumber] = index
            subformulas.append(formula)
            sumabundances.append([mass * prob, prob])
        if prob < cutoff:
            break

        tsv = tuple(v)
        if tsv not in vectorpool:
            ind = bisect.bisect(probabilityranking, r)
            probabilityranking.insert(ind, r)
            multinomialpath.insert(ind, (r, p, m, e, v))
            vectorpool.add(tsv)

        checkedcombos = set()
        for path in multinomialpath.copy():
            multielement = False
            match path[1]:
                case tuple():
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
                            sef = ''
                            for n, c in enumerate(sv):
                                if c > 0:
                                    sef += f'{fragmentpositions[se][n]}({c})'
                            seformulas.append(sef)
                            multipath.append((sr, sp, sm, se, sv))
                    checkformula = ''.join((sorted(seformulas)))
                    if checkformula in checkedcombos:
                        continue
                    else:
                        checkedcombos.add(checkformula)
                    if len(multipath) == 0:
                        continue
                case _:
                    sr, sp, sm, se, sv = path
                    sef = ''.join((f'{se}{str(n)}{(val)}' for n, val in enumerate(sv)))
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
                    seformula = ''
                    newprob = 1
                    newmass = 0
                    newmassnum = 0
                    newiter = {k: v for k, v in baseiter.items() if k not in sepool}
                    newiter[e] = [(r, p, m, e, v)]
                    for ir, ip, im, ie, iv in multipath:
                        newiter[ie] = [(ir, ip, im, ie, iv)]
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    iso = fragmentpositions[ie][n]
                                    newmassnum += int(iso[1:]) * c
                                    seformula += f'{iso}({c})'
                            newprob *= ip
                            newmass += im
                else:
                    newiter = {k: v for k, v in baseiter.items() if k != se}
                    newiter[se] = [(sr, sp, sm, se, sv)]
                    seformula = ''
                    newprob = 1
                    newmass = 0
                    newmassnum = 0
                    for b in sorted(newiter):
                        for ir, ip, im, ie, iv in newiter[b]:
                            for n, c in enumerate(iv):
                                if c > 0:
                                    iso = fragmentpositions[ie][n]
                                    newmassnum += int(iso[1:]) * c
                                    seformula += f'{iso}({c})'
                            newprob *= ip
                            newmass += im
                if newprob >= cutoff:
                    try:
                        index = massnumberindices[newmassnum]
                        subformulas[index] += '-' + seformula
                        sumabundances[index][0] += newmass * newprob
                        sumabundances[index][1] += newprob
                    except KeyError:
                        index = len(massnumberindices)
                        massnumberindices[newmassnum] = index
                        subformulas.append(seformula)
                        sumabundances.append([newmass * newprob, newprob])
                    if multielement:
                        ind = bisect.bisect(probabilityranking, newratio)
                        probabilityranking.insert(ind, newratio)
                        multinomialpath.insert(ind, (newratio, *multipath))
                    else:
                        newmulti = []
                        tsv = tuple(sv)
                        if tsv not in vectorpool:
                            newmulti.append((sr, sp, sm, se, sv))
                            vectorpool.add(tsv)
                        tvv = tuple(v)
                        if tvv not in vectorpool:
                            newmulti.append((r, p, m, e, v))
                            vectorpool.add(tvv)
                        if newmulti:
                            ind = bisect.bisect(probabilityranking, newratio)
                            probabilityranking.insert(ind, newratio)
                            multinomialpath.insert(ind, (newratio, *newmulti))
            else:
                break

    subformulas = np.array(subformulas, dtype='S')
    massesandabundances = np.array(sumabundances)
    massesandabundances[:, 0] /= massesandabundances[:, 1]
    order = massesandabundances[:, 1].argsort()[::-1]
    subformulas = subformulas[order].tolist()
    massesandabundances = massesandabundances[order]
    return subformulas, massesandabundances


def fragmentation_compositions(aminocomps, nfrags, cfrags, seq):
    """b/y (etc.) fragment elemental compositions for a peptide ``seq``."""
    fragments = {}

    fragcomp_n = {}
    for n, aa in enumerate(seq[:-1]):
        aa_composition = aminocomps[aa]
        for k in aa_composition:
            fragcomp_n[k] = fragcomp_n.get(k, 0) + aa_composition.get(k, 0)
        for ion, modcomp in nfrags.items():
            fragment_composition = fragcomp_n.copy()
            for k in modcomp:
                fc = fragment_composition.get(k, 0) + modcomp.get(k, 0)
                if fc > 0:
                    fragment_composition[k] = fc
                else:
                    fragment_composition.pop(k, None)
            fragments[ion + str(n + 1)] = fragment_composition

    fragcomp_c = {}
    for n, aa in enumerate(seq[::-1][:-1]):
        aa_composition = aminocomps[aa]
        for k in aa_composition:
            fragcomp_c[k] = fragcomp_c.get(k, 0) + aa_composition.get(k, 0)
        for ion, modcomp in cfrags.items():
            fragment_composition = fragcomp_c.copy()
            for k in modcomp:
                fc = fragment_composition.get(k, 0) + modcomp.get(k, 0)
                if fc > 0:
                    fragment_composition[k] = fc
                else:
                    fragment_composition.pop(k, None)
            fragments[ion + str(n + 1)] = fragment_composition

    return fragments


# ---- high-level entry point ------------------------------------------------

def _isoprobs_from_subformula(subformula):
    """Parse a precursor isotopologue subformula (``C12(40)C13(2)...``) into the
    per-element isotope-probability map the fragment walk expects (the reference
    main loop's ``isoprobs`` construction)."""
    isocounts = set()
    competing = set()
    competitors = {}
    isosums = {}
    for ss in subformula.split(')')[:-1]:
        iso, c = ss.split('(')
        c = int(c)
        e = iso[0]
        if e in isocounts:
            competing.add(e)
            competitors[e][iso] = c
            isosums[e] += c
        else:
            isocounts.add(e)
            competitors[e] = {iso: c}
            isosums[e] = c
    isoprobs = defaultdict(dict)
    for e, v in competitors.items():
        if e in competing:
            for iso, c in v.items():
                isoprobs[e][iso] = c / isosums[e]
        else:
            for iso in v:
                isoprobs[e][iso] = 1
    return isoprobs


def _fragment_main_masses(subformula, fragments, dividingthreshold):
    """Most-abundant fragment isotopomer neutral mass per ion, given a parent
    precursor isotopologue ``subformula`` (the reference main-loop body)."""
    isoprobs = _isoprobs_from_subformula(subformula)
    out = {}
    for ion, fragcomp in fragments.items():
        elementalorganizer = {}
        fragmentpositions = {}
        ok = True
        for e, c in fragcomp.items():
            if e in isoprobs:
                fragprobs = isoprobs[e]
                if len(fragprobs) > 1:
                    elementlist, positions = fragment_element_binomial_walk(
                        dividingthreshold, e, c, fragprobs)
                    elementalorganizer[e] = elementlist.copy()
                    fragmentpositions[e] = positions
                else:
                    iso = list(fragprobs)[0]
                    elementalorganizer[e] = [[-1, 1, elementalmasses[iso] * c, e, [c]]]
                    fragmentpositions[e] = {0: iso}
            else:
                # The parent isotopologue doesn't constrain this element (it had
                # no atoms of it) -> fall back to its monoisotope.
                mk = chem.MONOISOTOPIC_KEYS.get(e)
                if mk is None:
                    ok = False
                    break
                elementalorganizer[e] = [[-1, 1, elementalmasses[mk] * c, e, [c]]]
                fragmentpositions[e] = {0: mk}
        if not ok or not elementalorganizer:
            continue
        _formulas, ma = fragment_descending_partial_products(
            dividingthreshold, elementalorganizer, fragmentpositions)
        out[ion] = float(ma[0][0])   # most-abundant isotopomer neutral mass
    return out


def nearest_neighbors_ppm(baselist, flylist, ppm):
    """Port of ``nearest_neighbors_ppm_tolerance`` (peptidefragmentscoring.py).

    For each base m/z (theoretical, sorted ascending) finds the nearest fly m/z
    (experimental peaks, sorted ascending) within ``ppm``. Returns
    ``{base_index: [fly_index, ...]}`` (two indices only on an exact tie)."""
    base = np.asarray(baselist, dtype=float)
    fly = np.asarray(flylist, dtype=float)
    indices = {}
    if base.size == 0 or fly.size == 0:
        return indices
    ppmmod = ppm / 1e6
    rights = np.searchsorted(fly, base)
    for bn, rightfn in enumerate(rights.tolist()):
        b = base[bn]
        btol = b * ppmmod
        bmin = b - btol
        bmax = b + btol
        leftfn = rightfn - 1
        left = leftfn >= 0 and bmin < fly[leftfn] < bmax
        right = rightfn < fly.size and bmin < fly[rightfn] < bmax
        if left and right:
            leftdist = b - fly[leftfn]
            rightdist = fly[rightfn] - b
            if leftdist < rightdist:
                indices[bn] = [leftfn]
            elif rightdist < leftdist:
                indices[bn] = [rightfn]
            else:
                indices[bn] = [leftfn, rightfn]
        elif left:
            indices[bn] = [leftfn]
        elif right:
            indices[bn] = [rightfn]
    return indices


def annotate_spectrum(entries, peak_mz, peak_int, ppm=20.0):
    """Match theoretical fragment ions to the REAL MS2 peaks and split the
    spectrum into matched / unmatched actual peaks.

    Theoretical ``entries`` (from :func:`peptide_fragment_ions`) are matched to
    the experimental peaks via :func:`nearest_neighbors_ppm` at the search ppm.
    Returns ``(matched, (unmatched_mz, unmatched_int))`` where ``matched`` is a
    list of ``{mz, intensity, ion, charge, isotopes}`` at the EXPERIMENTAL peaks
    that were hit (theoretical ions with no match are simply dropped)."""
    fly = np.asarray(peak_mz, dtype=float)
    fint = np.asarray(peak_int, dtype=float)
    if fly.size == 0 or fint.size != fly.size:
        return [], (np.array([]), np.array([]))
    order = fly.argsort()
    fly = fly[order]
    fint = fint[order]

    ents = sorted(entries, key=lambda e: e["mz"])
    base = [e["mz"] for e in ents]
    nn = nearest_neighbors_ppm(base, fly, ppm)

    flyann = defaultdict(list)
    for bi, flist in nn.items():
        for fi in flist:
            flyann[fi].append(ents[bi])

    matched = []
    for fi in sorted(flyann):
        anns = flyann[fi]
        rep = min(anns, key=lambda a: abs(a["mz"] - fly[fi]))
        isos = sorted({a["isotope_index"] for a in anns if a["ion"] == rep["ion"]})
        matched.append({"mz": float(fly[fi]), "intensity": float(fint[fi]),
                        "ion": rep["ion"], "charge": rep["charge"], "isotopes": isos})

    matched_set = set(flyann)
    keep = np.array([i for i in range(fly.size) if i not in matched_set], dtype=int)
    if keep.size:
        return matched, (fly[keep], fint[keep])
    return matched, (np.array([]), np.array([]))


def match_fragment_ions(entries, peak_mz, peak_int, ppm=20.0):
    """Match theoretical fragment ions to an MS2 spectrum.

    Collapses the per-charge ``entries`` from :func:`peptide_fragment_ions` to
    one result per (ion, precursor-isotope), preferring the charge state that
    matched the most intense peak. Returns
    ``[{ion, isotope_index, charge, matched, mz, theo_mz, intensity}]`` —
    ``matched`` False means the ion is expected but absent (plot it red).
    """
    peaks = np.asarray(peak_mz, dtype=float)
    ints = np.asarray(peak_int, dtype=float)
    if peaks.size and ints.size == peaks.size:
        order = peaks.argsort()
        peaks = peaks[order]
        ints = ints[order]
    else:
        ints = None

    # Different precursor isotopes often yield the SAME fragment m/z (for a small
    # fragment the extra neutron usually lands in the complementary half), so key
    # by (ion, charge, rounded m/z) and collect the contributing isotope indices
    # -- one plotted line per distinct fragment, labelled with the isotopes it
    # could have descended from.
    groups = defaultdict(lambda: {"isos": set(), "mz": None, "charge": None,
                                  "ion": None})
    for e in entries:
        key = (e["ion"], e["charge"], round(e["mz"], 3))
        g = groups[key]
        g["ion"] = e["ion"]
        g["charge"] = e["charge"]
        g["mz"] = e["mz"]
        g["isos"].add(e["isotope_index"])

    out = []
    for g in groups.values():
        mz = g["mz"]
        tol = mz / 1e6 * ppm
        matched = False
        match_mz = mz
        intensity = 0.0
        if peaks.size:
            i = int(np.searchsorted(peaks, mz))
            for j in (i - 1, i):
                if 0 <= j < peaks.size and abs(peaks[j] - mz) <= tol:
                    inten = float(ints[j]) if ints is not None else 0.0
                    if not matched or inten > intensity:
                        matched = True
                        match_mz = float(peaks[j])
                        intensity = inten
        out.append({
            "ion": g["ion"],
            "charge": g["charge"],
            "isotopes": sorted(g["isos"]),
            "matched": matched,
            "mz": match_mz if matched else mz,
            "theo_mz": mz,
            "intensity": intensity,
        })
    return out


def peptide_fragment_ions(seq, charge, iso_low, iso_high, ions="by",
                          frag_charges=(1, 2),
                          dividingthreshold=DEFAULT_DIVIDING_THRESHOLD,
                          subisotopomericdepth=DEFAULT_SUBISO_DEPTH):
    """Theoretical b/y fragment ions for the precursor isotopes within the MS2
    isolation window ``[iso_low, iso_high]``.

    Returns ``[{ion, isotope_index, neutral_mass, mz}]``: ``isotope_index`` is
    the precursor isotope the fragment descended from (0 = monoisotopic), ``mz``
    is the singly/doubly-charged fragment m/z (one entry per ``frag_charges``).
    """
    seq = seq or ""
    if len(seq) < 2:
        return []
    charge = max(1, int(charge or 1))

    atomiccomposition = chem.peptide_atomic_composition(seq)
    try:
        subformulas, ma = distribution_generation(dividingthreshold, atomiccomposition)
    except Exception:
        return []
    subformulas = [s.decode() if isinstance(s, bytes) else s for s in subformulas]
    masses, abundances = ma

    # Group subformulas by nominal mass number (the precursor isotope group).
    groupmasses = defaultdict(list)
    groupabund = defaultdict(list)
    groupidx = defaultdict(list)
    for n, s in enumerate(subformulas):
        massnumber = 0
        for ss in s.split(")")[:-1]:
            iso, c = ss.split("(")
            massnumber += int(iso[1:]) * int(c)
        groupmasses[massnumber].append(masses[n])
        groupabund[massnumber].append(abundances[n])
        groupidx[massnumber].append(n)

    # Per group: mean mass + the subisotopomers within subisotopomericdepth.
    groups = []   # (massnumber, mean_mass, [global subformula indices])
    for mn in groupmasses:
        m = groupmasses[mn]
        a = groupabund[mn]
        gi = groupidx[mn]
        total = sum(a)
        weighted = 0.0
        cumulative = 0.0
        subinds = []
        cont = True
        for localn, (sm, sa) in sorted(enumerate(zip(m, a)), key=lambda x: -x[1][1]):
            weighted += sm * sa
            cumulative += sa
            cumpercent = cumulative / total if total else 1.0
            if cumpercent <= subisotopomericdepth:
                subinds.append(localn)
            elif cont and cumpercent >= subisotopomericdepth:
                subinds.append(localn)
                cont = False
        groups.append((mn, weighted / total if total else m[0],
                       [gi[li] for li in subinds]))

    # Isotope index = rank by mass (0 = monoisotopic). Select the groups whose
    # precursor m/z falls inside the MS2 isolation window.
    groups.sort(key=lambda g: g[1])
    fragments = fragmentation_compositions(
        chem.AMINO_ACID_COMPOSITION,
        {k: chem.N_FRAGMENT_COMPOSITIONS[k] for k in ions if k in chem.N_FRAGMENT_COMPOSITIONS},
        {k: chem.C_FRAGMENT_COMPOSITIONS[k] for k in ions if k in chem.C_FRAGMENT_COMPOSITIONS},
        seq)

    out = []
    seen = set()
    for isotope_index, (_mn, mean_mass, subidxs) in enumerate(groups):
        prec_mz = (mean_mass + proton * charge) / charge
        # Always generate the monoisotopic (index 0) fragments so a candidate is
        # never zeroed out just because its precursor isotope fell outside the
        # MS2 isolation window (e.g. a different-charge candidate); the window
        # only ADDS the higher precursor isotopes (M+1, M+2, ...) that were
        # co-isolated. This keeps Table 2 coverage and the panel-3 annotation
        # consistent for every candidate.
        if isotope_index != 0 and not (iso_low <= prec_mz <= iso_high):
            continue
        if not subidxs:
            continue
        # Use the most-abundant subisotopomer of this precursor isotope as the
        # representative parent for the fragment isotope walk.
        subformula = subformulas[subidxs[0]]
        try:
            main = _fragment_main_masses(subformula, fragments, dividingthreshold)
        except Exception:
            continue
        for ion, neutral in main.items():
            key = (ion, isotope_index)
            if key in seen:
                continue
            seen.add(key)
            for z in frag_charges:
                out.append({
                    "ion": ion,
                    "isotope_index": isotope_index,
                    "charge": z,
                    "neutral_mass": neutral,
                    "mz": (neutral + proton * z) / z,
                })
    return out
