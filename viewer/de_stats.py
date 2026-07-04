"""Differential-expression statistics for the Quantitative Comparisons tab.

Pure-Python (no scipy/pandas dependency — the runtime only guarantees numpy, and
these are small enough to keep dependency-free and unit-testable). Everything
operates on *log2* quantities: LFQ intensities are log-normal, so log2 makes the
fold change additive and the t-test's normality assumption reasonable.

Exposed:
    log2_values(quantities)          -> [log2(q) for q>0]
    welch_ttest(a, b)                -> (t, df, p_two_sided)
    paired_ttest(pairs)              -> (t, df, p_two_sided)
    benjamini_hochberg(pvalues)      -> [adjusted p, aligned with input]
    differential_expression(...)     -> per-feature DE records (see below)

The t-distribution survival function is computed from the regularized
incomplete beta (Numerical Recipes ``betai``/``betacf``), so no scipy is needed.
"""

import math


# --------------------------------------------------------------------------- #
# incomplete beta  ->  Student-t two-sided p-value
# --------------------------------------------------------------------------- #

def _betacf(a, b, x, itmax=200, eps=3.0e-12):
    """Continued-fraction expansion for the incomplete beta (Numerical Recipes)."""
    fpmin = 1.0e-30
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < fpmin:
        d = fpmin
    d = 1.0 / d
    h = d
    for m in range(1, itmax + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < fpmin:
            d = fpmin
        c = 1.0 + aa / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < eps:
            break
    return h


def _betai(a, b, x):
    """Regularized incomplete beta function I_x(a, b)."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    lbeta = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    bt = math.exp(lbeta)
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_sf_two_sided(t, df):
    """Two-sided p-value for a Student-t statistic with ``df`` degrees of freedom."""
    if df <= 0 or not math.isfinite(t):
        return float("nan")
    if t == 0.0:
        return 1.0
    x = df / (df + t * t)
    return _betai(df / 2.0, 0.5, x)


# --------------------------------------------------------------------------- #
# basic moments
# --------------------------------------------------------------------------- #

def _mean(xs):
    return sum(xs) / len(xs)


def _var(xs):
    """Sample variance (n-1). Returns 0.0 for length-1 input."""
    n = len(xs)
    if n < 2:
        return 0.0
    m = _mean(xs)
    return sum((x - m) ** 2 for x in xs) / (n - 1)


def log2_values(quantities):
    """log2 of every strictly-positive quantity (zeros/blanks are treated as
    missing and dropped, not floored to a pseudocount)."""
    out = []
    for q in quantities:
        if q is None:
            continue
        try:
            q = float(q)
        except (TypeError, ValueError):
            continue
        if q > 0.0:
            out.append(math.log2(q))
    return out


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #

def welch_ttest(a, b):
    """Welch's unequal-variance t-test. ``a``/``b`` are log2 value lists.

    Returns ``(t, df, p)``; ``p`` is NaN when either group has <2 usable values
    or both variances are zero."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return (float("nan"), float("nan"), float("nan"))
    va, vb = _var(a), _var(b)
    se2 = va / na + vb / nb
    if se2 <= 0.0:
        return (float("nan"), float("nan"), float("nan"))
    t = (_mean(a) - _mean(b)) / math.sqrt(se2)
    df_num = se2 * se2
    df_den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1)
    df = df_num / df_den if df_den > 0 else float("nan")
    return (t, df, t_sf_two_sided(t, df))


def paired_ttest(pairs):
    """Paired t-test over ``pairs`` = list of ``(a_log2, b_log2)`` for matched
    samples. Returns ``(t, df, p)``."""
    diffs = [a - b for a, b in pairs]
    n = len(diffs)
    if n < 2:
        return (float("nan"), float("nan"), float("nan"))
    sd = math.sqrt(_var(diffs))
    if sd <= 0.0:
        return (float("nan"), float("nan"), float("nan"))
    t = _mean(diffs) / (sd / math.sqrt(n))
    df = n - 1
    return (t, df, t_sf_two_sided(t, df))


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg FDR adjustment. NaN p-values pass through as NaN and
    are excluded from the ranking. Returns a list aligned with ``pvalues``."""
    indexed = [(i, p) for i, p in enumerate(pvalues)
               if p is not None and not math.isnan(p)]
    m = len(indexed)
    adjusted = [float("nan")] * len(pvalues)
    if m == 0:
        return adjusted
    indexed.sort(key=lambda ip: ip[1])
    prev = 1.0
    # Walk from largest p to smallest, enforcing monotonicity.
    for rank in range(m, 0, -1):
        i, p = indexed[rank - 1]
        val = min(prev, p * m / rank)
        adjusted[i] = val
        prev = val
    return adjusted


# --------------------------------------------------------------------------- #
# feature-level DE
# --------------------------------------------------------------------------- #

def differential_expression(features, samples_a, samples_b, matrix,
                            paired_pairs=None, min_replicates=2):
    """Run DE for many features at once.

    ``features``      : iterable of feature keys (peptide/protein ids)
    ``samples_a/b``   : filenames belonging to group A / group B
    ``matrix``        : {feature: {filename: quantity}}
    ``paired_pairs``  : optional list of ``(a_filename, b_filename)`` tuples; when
                        given a paired t-test is used over these matched samples.
    ``min_replicates``: minimum usable (positive) values required per group.

    Returns a list of dicts (one per feature that met ``min_replicates`` in both
    groups), each with: feature, mean_a, mean_b, log2fc, n_a, n_b, t, df, p,
    and (after the fact) fdr — sorted by ascending p (NaN p last). log2fc is
    ``mean_a - mean_b`` in log2 space, i.e. positive means higher in A.
    """
    records = []
    for feat in features:
        row = matrix.get(feat, {})
        if paired_pairs:
            pairs = []
            for fa, fb in paired_pairs:
                qa, qb = row.get(fa), row.get(fb)
                la = log2_values([qa])
                lb = log2_values([qb])
                if la and lb:
                    pairs.append((la[0], lb[0]))
            a_vals = [p[0] for p in pairs]
            b_vals = [p[1] for p in pairs]
            if len(pairs) < min_replicates:
                continue
            t, df, p = paired_ttest(pairs)
        else:
            a_vals = log2_values(row.get(f) for f in samples_a)
            b_vals = log2_values(row.get(f) for f in samples_b)
            if len(a_vals) < min_replicates or len(b_vals) < min_replicates:
                continue
            t, df, p = welch_ttest(a_vals, b_vals)

        mean_a = _mean(a_vals) if a_vals else float("nan")
        mean_b = _mean(b_vals) if b_vals else float("nan")
        records.append({
            "feature": feat,
            "mean_a": mean_a,
            "mean_b": mean_b,
            "log2fc": mean_a - mean_b,
            "n_a": len(a_vals),
            "n_b": len(b_vals),
            "t": t,
            "df": df,
            "p": p,
        })

    fdr = benjamini_hochberg([r["p"] for r in records])
    for r, q in zip(records, fdr):
        r["fdr"] = q

    def sort_key(r):
        p = r["p"]
        return (1, 0.0) if (p is None or math.isnan(p)) else (0, p)

    records.sort(key=sort_key)
    return records
