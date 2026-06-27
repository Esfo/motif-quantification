"""Envelope-first isotope-distribution / charge assignment.

This replaces the edge-first charge call (pair spacing -> rounded charge ->
per-charge path walking -> length-weighted greedy claiming). That structure
decided charge from a single local adjacent spacing before the whole isotope
cluster was ever fitted, which let higher charges win by forming longer local
paths in a dense feature map.

Here charge is decided at the *envelope* level:

  1. For each seed feature, derive the set of plausible charges from the spacing
     to its coeluting m/z neighbours (allowing isotope-index gaps q>1, so a
     missing M+1 is not misread as a different charge). No hard min/max charge.

  2. For each candidate charge, fit the whole isotope lattice
         mz(k) = mono_mz + k * C13 / z
     by walking forward and backward from the seed, occupying observed features
     at expected positions and tolerating a couple of missing rungs.

  3. Fit the monoisotope offset (which observed rung is M+0) by averagine fit,
     instead of assuming the leftmost member is M+0.

  4. Score the *whole* envelope — m/z lattice residual, averagine intensity
     agreement, chromatographic trace-shape consensus, and missing-expected-peak
     penalty — with NO raw path-length reward. A clean 3-peak 2+ envelope is no
     longer beaten by a longer high-charge path just because it has more rungs.

  5. Let complete explanations compete for features (a feature belongs to one
     envelope); the winner is the charge whose envelope actually fits best.

Output rows are the same dict shape the sqlite writer expects, plus per-envelope
and per-member evidence fields (used by the GUI / schema in a later phase).
"""
import math

import numpy as np

PROTON = 1.007276554940804
C13_DELTA = 1.00335483507
AVERAGINE_LAMBDA_PER_DA = 0.000594  # Poisson mean of the averagine isotope model

MAX_ISOTOPE_GAP = 2     # consecutive missing rungs tolerated while growing
MAX_OFFSET_SEARCH = 3   # monoisotope offsets (M+0..M+3) considered for the first rung
MAX_Q = 2               # isotope-index gap considered when deriving charge from a pair


def averagine_envelope(neutral_mass, n_peaks):
    """Normalised (max=1) averagine isotope intensities for indices 0..n_peaks-1."""
    lam = max(neutral_mass, 1.0) * AVERAGINE_LAMBDA_PER_DA
    ks = np.arange(n_peaks)
    logp = -lam + ks * math.log(lam) - np.array([math.lgamma(k + 1) for k in ks])
    p = np.exp(logp)
    m = p.max()
    return p / m if m > 0 else p


class EnvelopeContext:
    """Vectorised feature arrays + trace lookup shared by all candidate fits."""

    def __init__(self, features, traces, config_dict):
        self.cfg = config_dict
        self.mz = np.asarray([f.mz_mean for f in features], dtype=np.float64)
        self.rt_start = np.asarray([f.rt_start for f in features], dtype=np.float64)
        self.rt_apex = np.asarray([f.rt_apex for f in features], dtype=np.float64)
        self.rt_end = np.asarray([f.rt_end for f in features], dtype=np.float64)
        self.ms1_start = np.asarray([f.ms1_start for f in features], dtype=np.int64)
        self.ms1_apex = np.asarray([f.ms1_apex for f in features], dtype=np.int64)
        self.ms1_end = np.asarray([f.ms1_end for f in features], dtype=np.int64)
        self.height = np.asarray([f.height for f in features], dtype=np.float64)
        self.area = np.asarray([f.area for f in features], dtype=np.float64)
        self.n_points = np.asarray([f.n_points for f in features], dtype=np.int64)
        self.traces = traces  # list indexed by feature_id: {scans,rts,mzs,intensities}
        self.order = np.argsort(self.mz)
        self.sorted_mz = self.mz[self.order]
        self._sim_cache = {}

    # ---- m/z neighbour search ------------------------------------------------

    def mz_tol(self, mz_value):
        return max(
            self.cfg["isotope_mz_abs"],
            mz_value * self.cfg["isotope_mz_ppm"] / 1_000_000.0,
        )

    def features_near(self, target_mz):
        """Feature ids whose m/z is within tolerance of target_mz."""
        tol = self.mz_tol(target_mz)
        lo = np.searchsorted(self.sorted_mz, target_mz - tol, side="left")
        hi = np.searchsorted(self.sorted_mz, target_mz + tol, side="right")
        return [int(self.order[k]) for k in range(lo, hi)]

    def features_in_range(self, lo_mz, hi_mz):
        lo = np.searchsorted(self.sorted_mz, lo_mz, side="left")
        hi = np.searchsorted(self.sorted_mz, hi_mz, side="right")
        return [int(self.order[k]) for k in range(lo, hi)]

    def neighbours_above(self, seed, max_dmz):
        seed_mz = self.mz[seed]
        lo = np.searchsorted(self.sorted_mz, seed_mz, side="right")
        hi = np.searchsorted(self.sorted_mz, seed_mz + max_dmz, side="right")
        out = []
        for k in range(lo, hi):
            j = int(self.order[k])
            if j != seed:
                out.append(j)
        return out

    # ---- chromatographic trace similarity -----------------------------------

    def coelutes(self, a, b):
        """RT-window overlap gate (cheap pre-filter before trace cosine)."""
        overlap = min(self.rt_end[a], self.rt_end[b]) - max(self.rt_start[a], self.rt_start[b])
        return overlap > 0

    def rt_overlap_ratio(self, a, b):
        overlap = min(self.rt_end[a], self.rt_end[b]) - max(self.rt_start[a], self.rt_start[b])
        union = max(self.rt_end[a], self.rt_end[b]) - min(self.rt_start[a], self.rt_start[b])
        return max(0.0, float(overlap / union)) if union > 0 else 0.0

    def coelution_score(self, a, b):
        """Membership score for attaching feature b to a's lattice.

        Normally the trace-shape cosine. When short_trace_fallback is on (the
        recovery pass) and either feature is too short for a reliable cosine, fall
        back to apex-RT proximity + window overlap and take the better of the two —
        so a weak/short isotope line that genuinely coelutes is not discarded just
        because a handful of noisy scans gave it a low cosine.
        """
        sim = self.trace_similarity(a, b)
        if not self.cfg.get("short_trace_fallback", False):
            return sim
        if min(int(self.n_points[a]), int(self.n_points[b])) >= self.cfg.get("short_trace_len", 5):
            return sim
        apex_gap = abs(self.rt_apex[a] - self.rt_apex[b])
        width = max(self.rt_end[a] - self.rt_start[a], self.rt_end[b] - self.rt_start[b], 1e-9)
        prox = max(0.0, 1.0 - apex_gap / width)
        proximity = 0.5 * prox + 0.5 * self.rt_overlap_ratio(a, b)
        return max(sim, proximity)

    def trace_similarity(self, a, b):
        """Cosine of the two elution traces aligned on shared scans (0..1).

        Real isotope members of one analyte share a chromatographic shape; an
        accidental coeluting peak usually does not. Missing scans count as zero
        intensity so partial overlap is penalised.
        """
        if a == b:
            return 1.0
        key = (a, b) if a < b else (b, a)
        cached = self._sim_cache.get(key)
        if cached is not None:
            return cached
        ta = self.traces[a]
        tb = self.traces[b]
        sa, ia = ta["scans"], ta["intensities"]
        sb, ib = tb["scans"], tb["intensities"]
        lo = int(min(sa.min(), sb.min()))
        hi = int(max(sa.max(), sb.max()))
        n = hi - lo + 1
        if n <= 0:
            self._sim_cache[key] = 0.0
            return 0.0
        va = np.zeros(n, dtype=np.float64)
        vb = np.zeros(n, dtype=np.float64)
        va[sa - lo] = ia
        vb[sb - lo] = ib
        na = np.linalg.norm(va)
        nb = np.linalg.norm(vb)
        sim = float(np.dot(va, vb) / (na * nb)) if na > 0 and nb > 0 else 0.0
        sim = max(0.0, min(1.0, sim))
        self._sim_cache[key] = sim
        return sim


    def consensus_trace_scores(self, fids):
        """Build a consensus elution trace and score each member against it.

        A real isotope envelope shares one chromatographic shape, so each member
        should match the intensity-summed consensus. Returns (scores_per_member,
        mean_score). This is stronger than matching every member to an arbitrary
        seed, which could itself be a tail peak.
        """
        traces = [self.traces[f] for f in fids]
        lo = int(min(t["scans"].min() for t in traces))
        hi = int(max(t["scans"].max() for t in traces))
        n = hi - lo + 1
        if n <= 0:
            return [1.0] * len(fids), 1.0
        vecs = []
        consensus = np.zeros(n, dtype=np.float64)
        for t in traces:
            v = np.zeros(n, dtype=np.float64)
            v[t["scans"] - lo] = t["intensities"]
            vecs.append(v)
            consensus += v
        cn = np.linalg.norm(consensus)
        if cn <= 0:
            return [1.0] * len(fids), 1.0
        scores = []
        for v in vecs:
            vn = np.linalg.norm(v)
            scores.append(float(np.dot(v, consensus) / (vn * cn)) if vn > 0 else 0.0)
        scores = [max(0.0, min(1.0, s)) for s in scores]
        return scores, float(np.mean(scores))


def derive_charges(ctx, seed):
    """Plausible charges for an envelope seeded at `seed`.

    Read off the spacing to each coeluting neighbour just above the seed, for
    isotope-index gaps q in 1..MAX_Q, so a missing M+1 (q=2 spacing) is not
    misread as a different charge. Returns a set of candidate integer charges.
    """
    charges = set()
    max_dmz = MAX_Q * C13_DELTA + ctx.mz_tol(ctx.mz[seed] + C13_DELTA)
    for j in ctx.neighbours_above(seed, max_dmz):
        dmz = ctx.mz[j] - ctx.mz[seed]
        if dmz <= 0:
            continue
        if not ctx.coelutes(seed, j):
            continue
        for q in range(1, MAX_Q + 1):
            z = int(round(q * C13_DELTA / dmz))
            if z < 1:
                continue
            expected = q * C13_DELTA / z
            if abs(dmz - expected) <= ctx.mz_tol(ctx.mz[j]):
                charges.add(z)
    return charges


def grow_lattice(ctx, seed, z):
    """Occupy observed features on the z-lattice anchored at the seed.

    Returns {relative_index: feature_id} with the seed at index 0. Walks forward
    and backward, picking at each expected rung the coeluting feature with the
    best trace similarity to the seed, tolerating up to MAX_ISOTOPE_GAP misses.
    """
    seed_mz = ctx.mz[seed]
    occupied = {0: seed}
    spacing = C13_DELTA / z
    min_sim = ctx.cfg.get("min_trace_similarity", 0.5)

    for direction in (1, -1):
        misses = 0
        k = direction
        while misses <= MAX_ISOTOPE_GAP:
            target = seed_mz + k * spacing
            if target <= 0:
                break
            best_j, best_sim = None, min_sim
            for j in ctx.features_near(target):
                if j == seed or j in occupied.values():
                    continue
                if not ctx.coelutes(seed, j):
                    continue
                sim = ctx.coelution_score(seed, j)
                if sim >= best_sim:
                    best_sim = sim
                    best_j = j
            if best_j is None:
                misses += 1
            else:
                occupied[k] = best_j
                misses = 0
            k += direction
    return occupied


def score_envelope(ctx, occupied, z):
    """Fit monoisotope offset and score the whole envelope.

    Returns the best (score_dict, members, mono_index) over the offset search, or
    None if the envelope is below the per-charge member floor.
    """
    rel_indices = sorted(occupied)
    fids = [occupied[r] for r in rel_indices]
    n_obs = len(fids)
    if n_obs < min_members_for_charge(z, ctx.cfg):
        return None

    mz = ctx.mz[fids]
    heights = ctx.height[fids]
    spacing = C13_DELTA / z
    base_rel = rel_indices[0]

    best = None
    for offset in range(MAX_OFFSET_SEARCH + 1):
        # The first observed rung is M+offset; isotope index of each member:
        iso_idx = np.array([r - base_rel + offset for r in rel_indices], dtype=np.int64)
        mono_mz = mz[0] - iso_idx[0] * spacing
        neutral_mass = (mono_mz - PROTON) * z
        if neutral_mass <= 0 or neutral_mass > ctx.cfg["max_neutral_mass"]:
            continue

        span = int(iso_idx.max()) + 1
        expected = averagine_envelope(neutral_mass, max(span, n_obs))

        # ---- m/z lattice residual ----
        lattice_mz = mono_mz + iso_idx * spacing
        resid = np.abs(mz - lattice_mz)
        tol = np.array([ctx.mz_tol(m) for m in mz])
        mz_score = float(np.mean(np.clip(1.0 - resid / tol, 0.0, 1.0)))

        # ---- averagine intensity agreement ----
        exp_obs = expected[iso_idx]
        obs_norm = heights / heights.max() if heights.max() > 0 else heights
        # scale expected to observed by least squares then cosine-like agreement
        denom = float(np.dot(exp_obs, exp_obs))
        scale = float(np.dot(obs_norm, exp_obs) / denom) if denom > 0 else 0.0
        pred = scale * exp_obs
        num = float(np.dot(obs_norm, pred))
        na = float(np.linalg.norm(obs_norm))
        nb = float(np.linalg.norm(pred))
        iso_score = max(0.0, num / (na * nb)) if na > 0 and nb > 0 else 0.0

        # ---- missing-expected-peak penalty (detection-aware) ----
        # An expected rung is only "missing" if it should have been visible:
        # scale the averagine to the observed members and compare against the
        # dimmest observed member (a proxy for the local detection floor). A rung
        # predicted below that floor is below detection and is not penalised.
        present = set(int(i) for i in iso_idx)
        obs_floor = float(obs_norm.min())
        pred_full = scale * expected  # expected scaled into observed-norm units
        missing_interior = 0
        for k in range(int(iso_idx.min()), int(iso_idx.max()) + 1):
            if k not in present and pred_full[k] > obs_floor:
                missing_interior += 1
        # monoisotope expectation: for low mass M+0 should be clearly visible
        if 0 not in present and pred_full[0] > max(obs_floor, 0.4):
            missing_interior += 1
        missing_score = 1.0 / (1.0 + missing_interior)

        # ---- trace-shape consensus ----
        member_trace_scores, trace_score = ctx.consensus_trace_scores(fids)

        # ---- interloper penalty (kills decimated aliases) ----
        # A coeluting feature inside the envelope's m/z span that does NOT sit on
        # the lattice is unexplained. The decisive case: the M+1/M+3 rungs of a
        # real z=2 envelope are interlopers under the z=1 (double-spacing) alias,
        # whereas the true z=2 fit has none. Penalising them makes the higher,
        # denser-but-consistent charge win without any charge prior.
        fid_set = set(int(f) for f in fids)
        lattice_mzs = lattice_mz
        interlopers = 0
        margin = spacing * 0.5
        for j in ctx.features_in_range(float(mz.min()) - 1e-6, float(mz.max()) + 1e-6):
            if j in fid_set:
                continue
            if not ctx.coelutes(fids[0], j):
                continue
            if ctx.coelution_score(fids[0], j) < ctx.cfg.get("min_trace_similarity", 0.5):
                continue
            jmz = ctx.mz[j]
            on_lattice = bool(np.any(np.abs(lattice_mzs - jmz) <= ctx.mz_tol(jmz)))
            # only count interlopers that fall strictly between rungs (not a tail
            # member we simply didn't grow into)
            if not on_lattice and float(np.min(np.abs(lattice_mzs - jmz))) < margin:
                interlopers += 1
        interloper_score = 1.0 / (1.0 + interlopers)

        total = (
            0.22 * mz_score
            + 0.34 * iso_score
            + 0.18 * trace_score
            + 0.13 * missing_score
            + 0.13 * interloper_score
        )

        if best is None or total > best[0]["total"]:
            members = [
                {
                    "feature_id": int(fids[i]),
                    "isotope_index": int(iso_idx[i]),
                    "member_score": float(min(1.0, max(0.0, 1.0 - resid[i] / tol[i]))),
                    "mz_residual": float(resid[i]),
                    "intensity_observed": float(obs_norm[i]),
                    "intensity_expected": float(exp_obs[i]),
                    "trace_score": float(member_trace_scores[i]),
                }
                for i in range(n_obs)
            ]
            best = (
                {
                    "total": float(total),
                    "mz_score": mz_score,
                    "iso_score": float(iso_score),
                    "trace_score": trace_score,
                    "missing_score": missing_score,
                    "interloper_score": float(interloper_score),
                    "n_interlopers": int(interlopers),
                    "neutral_mass": float(neutral_mass),
                    "mono_mz": float(mono_mz),
                    "mono_offset": int(offset),
                    "n_missing_interior": int(missing_interior),
                },
                members,
            )

    return best


def min_members_for_charge(charge, cfg):
    if charge == 1:
        return cfg["min_members_charge_one"]
    return cfg["min_distribution_members"]


def build_distributions_two_pass(features, traces, config_dict, progress=False):
    """Primary envelope pass + a relaxed recovery pass on the leftovers.

    The primary pass is unchanged (the clean, high-confidence set). The features
    no distribution claimed are then re-run through the same builder with relaxed
    gates (lower trace-similarity and envelope-score floors) and the short-trace
    coelution fallback, recovering weaker/shorter envelopes that are visually
    obvious in m/z-vs-RT but too faint for the strict pass. Recovered rows carry
    status='recovered' so they stay a separate, opt-in confidence tier.
    """
    primary = build_distributions_envelope(features, traces, config_dict, progress=progress)

    if not config_dict.get("enable_recovery", True):
        return primary

    claimed = {m["feature_id"] for row in primary for m in row["members"]}
    leftover = [i for i in range(len(features)) if i not in claimed]
    if not leftover:
        return primary

    sub_features = [features[i] for i in leftover]
    sub_traces = [traces[i] for i in leftover]

    relaxed = dict(config_dict)
    relaxed["min_trace_similarity"] = config_dict.get("recover_min_trace_similarity", 0.3)
    relaxed["min_envelope_score"] = config_dict.get("recover_min_envelope_score", 0.30)
    relaxed["short_trace_fallback"] = True

    rec_rows = build_distributions_envelope(sub_features, sub_traces, relaxed, progress=progress)
    for row in rec_rows:
        row["status"] = "recovered"
        for m in row["members"]:
            m["feature_id"] = leftover[m["feature_id"]]  # sub-index -> original feature_id

    if progress:
        print(f"recovery leftover_features={len(leftover)} recovered={len(rec_rows)}",
              file=__import__("sys").stderr, flush=True)

    return primary + rec_rows


def build_distributions_envelope(features, traces, config_dict, progress=False):
    """Envelope-first distribution builder. Returns a list of row dicts."""
    if not features:
        return []

    ctx = EnvelopeContext(features, traces, config_dict)
    n = len(features)
    min_total = config_dict.get("min_envelope_score", 0.45)

    candidates = []
    seen = set()
    for seed in range(n):
        for z in derive_charges(ctx, seed):
            occupied = grow_lattice(ctx, seed, z)
            if len(occupied) < min_members_for_charge(z, config_dict):
                continue
            fid_key = (z, frozenset(occupied.values()))
            if fid_key in seen:
                continue
            seen.add(fid_key)
            scored = score_envelope(ctx, occupied, z)
            if scored is None:
                continue
            ev, members = scored
            if ev["total"] < min_total:
                continue
            candidates.append((ev, members, z))

        if progress and seed > 0 and seed % 20000 == 0:
            print(f"envelope seeds={seed}/{n} candidates={len(candidates)}",
                  file=__import__("sys").stderr, flush=True)

    return _compete(ctx, candidates)


def _compete(ctx, candidates):
    """Local model competition: a feature belongs to one envelope.

    Candidates are taken best-total first; a candidate is kept only if none of
    its features are already claimed. Because `total` is charge-neutral (no
    length reward), the winner for a shared feature set is the charge whose
    envelope genuinely fits best — not the one that formed the longest path.

    A dropped candidate of a *different charge* that conflicts with a winner is
    recorded as that winner's runner-up: when the two scores are within
    ambiguity_margin the winner is marked 'ambiguous' instead of 'validated', so
    a near-tie between charges is reported as ambiguity rather than a forced call.
    """
    candidates.sort(key=lambda c: c[0]["total"], reverse=True)
    margin = ctx.cfg.get("ambiguity_margin", 0.05)
    winner_of = {}        # feature_id -> index into `kept`
    runner_up = {}        # kept index -> best conflicting different-charge total
    claimed = set()
    kept = []
    for ev, members, z in candidates:
        fids = [m["feature_id"] for m in members]
        conflict = [f for f in fids if f in claimed]
        if conflict:
            # record this as a runner-up for any winner it overlaps at a
            # different charge (the ambiguity signal)
            for f in conflict:
                wi = winner_of.get(f)
                if wi is not None and kept[wi][2] != z:
                    runner_up[wi] = max(runner_up.get(wi, 0.0), ev["total"])
            continue
        idx = len(kept)
        kept.append((ev, members, z))
        for f in fids:
            claimed.add(f)
            winner_of[f] = idx

    rows = []
    for idx, (ev, members, z) in enumerate(kept):
        fids = [m["feature_id"] for m in members]
        runner = runner_up.get(idx, 0.0)
        ambiguity_score = float(runner / ev["total"]) if ev["total"] > 0 else 0.0
        status = "ambiguous" if (ev["total"] - runner) < margin * ev["total"] and runner > 0 else "validated"

        path_array = np.asarray(fids, dtype=np.int64)
        apex_fid = int(path_array[np.argmax(ctx.height[path_array])])
        rows.append(
            {
                "charge": int(z),
                "neutral_mass": ev["neutral_mass"],
                "mono_mz": ev["mono_mz"],
                "rt_start": float(np.min(ctx.rt_start[path_array])),
                "rt_apex": float(ctx.rt_apex[apex_fid]),
                "rt_end": float(np.max(ctx.rt_end[path_array])),
                "ms1_start": int(np.min(ctx.ms1_start[path_array])),
                "ms1_apex": int(ctx.ms1_apex[apex_fid]),
                "ms1_end": int(np.max(ctx.ms1_end[path_array])),
                "n_members": len(fids),
                "score": ev["total"],
                "quality": ev["total"],
                "mz_score": ev["mz_score"],
                "iso_score": ev["iso_score"],
                "trace_score": ev["trace_score"],
                "missing_score": ev["missing_score"],
                "interloper_score": ev["interloper_score"],
                "mono_offset": ev["mono_offset"],
                "n_missing_interior": ev["n_missing_interior"],
                "n_interlopers": ev["n_interlopers"],
                "ambiguity_score": ambiguity_score,
                "status": status,
                "members": members,
            }
        )
    return rows
