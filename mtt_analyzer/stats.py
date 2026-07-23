"""Statistical functions for MTT analysis.

- mean / sd / sem / cv: basic descriptive statistics
- welch_t_test: Welchs t-test (no equal-variance assumption)
- holm_bonferroni: Holm-Bonferroni multiple-testing correction
"""

from __future__ import annotations

import math
from typing import Sequence, List, Dict, Tuple

import numpy as np


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


def _log_gamma(x):
    """Lanczos approximation of log(Gamma(x))."""
    c = [
        76.18009172947146,
        -86.50532032941677,
        24.01409824083091,
        -1.231739572450155,
        0.1208650973866179e-2,
        -0.5395239384953e-5,
    ]
    y = x
    tmp = x + 5.5
    tmp -= (x + 0.5) * math.log(tmp)
    ser = 1.000000000190015
    for j in range(6):
        y += 1
        ser += c[j] / y
    return -tmp + math.log(2.5066282746310005 * ser / x)


def incomplete_beta(a, b, x):
    """Regularized incomplete beta function I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    ln_beta = _log_gamma(a) + _log_gamma(b) - _log_gamma(a + b)
    front = math.exp(a * math.log(x) + b * math.log(1 - x) - ln_beta) / a
    f = 1.0
    cc = 1.0
    d = 0.0
    for i in range(200):
        m = i / 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1 + num * d
        if abs(d) < 1e-30:
            d = 1e-30
        cc = 1 + num / cc
        if abs(cc) < 1e-30:
            cc = 1e-30
        delta = cc / d
        f *= delta
        if abs(delta - 1) < 1e-10:
            break
    return front * (f - 1)


def welch_t_test(a, b):
    """Welch t-test (does not assume equal variance). Returns dict with t, df, p."""
    a1 = np.asarray(a, dtype=float)
    a2 = np.asarray(b, dtype=float)
    m1, m2 = float(np.mean(a1)), float(np.mean(a2))
    v1 = float(np.var(a1, ddof=1)) / len(a1) if len(a1) > 1 else 0.0
    v2 = float(np.var(a2, ddof=1)) / len(a2) if len(a2) > 1 else 0.0
    if v1 + v2 == 0:
        return {"t": 0.0, "df": 1.0, "p": 1.0}
    se = math.sqrt(v1 + v2)
    t = (m1 - m2) / se
    d_num = (v1 ** 2) / max(len(a1) - 1, 1) + (v2 ** 2) / max(len(a2) - 1, 1)
    df = 1.0 if d_num == 0 else (v1 + v2) ** 2 / d_num
    xx = df / (df + t * t)
    p = incomplete_beta(df / 2, 0.5, xx)
    return {"t": t, "df": df, "p": p}


def holm_bonferroni(p_values):
    """Holm-Bonferroni step-down correction."""
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    adjusted = [0.0] * n
    cum_min = 1.0
    for k, (orig_idx, p) in enumerate(indexed):
        adj = min(1.0, p * (n - k))
        cum_min = min(cum_min, adj)
        adjusted[orig_idx] = cum_min
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
