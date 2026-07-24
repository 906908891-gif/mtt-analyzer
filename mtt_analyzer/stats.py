"""Statistical functions for MTT analysis.

- mean / sd / sem / cv: basic descriptive statistics
- welch_t_test: Welch t-test (no equal-variance assumption)
- holm_bonferroni: Holm-Bonferroni multiple-testing correction
"""

from __future__ import annotations

import math
from typing import Sequence, List, Dict

import numpy as np


try:
    from scipy.special import betainc
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


def mean(values):
    arr = np.asarray(values, dtype=float)
    if len(arr) == 0:
        return float("nan")
    return float(np.mean(arr))


def sd(values):
    arr = np.asarray(values, dtype=float)
    if len(arr) < 2:
        return 0.0
    return float(np.std(arr, ddof=1))


def sem(values):
    arr = np.asarray(values, dtype=float)
    if len(arr) < 2:
        return 0.0
    return float(np.std(arr, ddof=1) / math.sqrt(len(arr)))


def cv(values):
    """Coefficient of variation (%)."""
    m = mean(values)
    if m == 0 or not math.isfinite(m):
        return float("nan")
    return sd(values) / abs(m) * 100.0


# ----- Incomplete beta function -----
# Standard implementation from Numerical Recipes (Press et al., 1992)
# After Lentz continued fraction expansion.

def incomplete_beta(a, b, x):
    """Regularized incomplete beta function I_x(a, b).

    Uses scipy.special.betainc when available; otherwise uses the NR
    continued-fraction implementation.

    Numerical verification: I_0.5(1, 1) = 0.5
    """
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    if _HAS_SCIPY:
        return float(betainc(a, b, x))
    # Numerical Recipes continued fraction (modified Lentz)
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def _betacf(a, b, x, max_iter=200):
    """Continued fraction for incomplete beta (NR Press et al.)."""
    EPS = 1e-12
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        # Even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0 / d
        h *= d * c
        # Odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            break
    return h


# ----- Welch t-test -----

def welch_t_test(a, b):
    """Welch t-test (does not assume equal variance).

    Returns a dict with keys: t, df, p.
    """
    a1 = np.asarray(a, dtype=float)
    a2 = np.asarray(b, dtype=float)
    if len(a1) < 2 or len(a2) < 2:
        return {"t": 0.0, "df": 1.0, "p": 1.0}
    m1, m2 = float(np.mean(a1)), float(np.mean(a2))
    v1 = float(np.var(a1, ddof=1)) / len(a1)
    v2 = float(np.var(a2, ddof=1)) / len(a2)
    if v1 + v2 == 0:
        return {"t": 0.0, "df": 1.0, "p": 1.0}
    se = math.sqrt(v1 + v2)
    t_val = (m1 - m2) / se
    d_num = (v1 ** 2) / (len(a1) - 1) + (v2 ** 2) / (len(a2) - 1)
    df = 1.0 if d_num == 0 else (v1 + v2) ** 2 / d_num
    xx = df / (df + t_val * t_val)
    p = incomplete_beta(df / 2, 0.5, xx)
    p = max(0.0, min(1.0, p))
    return {"t": t_val, "df": df, "p": p}


# ----- Holm-Bonferroni -----

def holm_bonferroni(p_values):
    """Holm-Bonferroni step-down correction.

    Correct algorithm (step-down with running MAX):
      1. Sort p-values ascending
      2. adj_k = p_(k) * (n - k + 1)
      3. Enforce monotonicity with running MAX (not min)

    Reference: Holm, S. (1979).
    """
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    cum_max = 0.0
    for k, (orig_idx, p) in enumerate(indexed):
        adj = min(1.0, float(p) * (n - k))
        cum_max = max(cum_max, adj)
        adjusted[orig_idx] = cum_max
    return adjusted


def significance_stars(p):
    """Return significance stars for a p-value."""
    if p is None or not math.isfinite(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    if p < 0.1:
        return "."
    return "ns"
