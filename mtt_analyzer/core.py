"""4PL (four-parameter logistic) dose-response curve fitting.

Model
-----
y = B + (T - B) / (1 + 10 ** ((x - logEC50) * h))

where:
    T  = Top (response at zero dose)
    B  = Bottom (response at infinite dose)
    h  = HillSlope (steepness)
    logEC50 = log10 of EC50 (concentration giving 50% effect)

This implementation uses scipy.optimize.curve_fit (Levenberg-Marquardt)
when available, and falls back to a pure-numpy implementation otherwise.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, asdict
from typing import Dict, Any

import numpy as np


try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


@dataclass
class FourPLResult:
    Top: float
    Bottom: float
    logEC50: float
    HillSlope: float
    EC50: float
    R2: float
    ssRes: float
    ssTot: float
    nPoints: int
    converged: bool
    seTop: float = float("nan")
    seBottom: float = float("nan")
    seLogEC50: float = float("nan")
    seHillSlope: float = float("nan")
    ci95Low: float = float("nan")
    ci95High: float = float("nan")
    iterations: int = 0
    method: str = "unknown"

    def to_dict(self):
        return asdict(self)


def _four_pl(x, T, B, logEC50, h):
    x = np.asarray(x, dtype=float)
    return B + (T - B) / (1.0 + np.power(10.0, (x - logEC50) * h))


def _initial_guess(log_x, y):
    sorted_lx = np.sort(log_x)
    n = len(sorted_lx)
    return [float(np.max(y)), float(np.min(y)), float(sorted_lx[n // 2]), 1.0]


def _fit_scipy(log_x, y):
    if len(log_x) < 4:
        return None
    p0 = _initial_guess(log_x, y)
    try:
        popt, pcov = curve_fit(
            _four_pl, log_x, y, p0=p0,
            bounds=([-float("inf"), -float("inf"), -20, 1e-4],
                    [float("inf"), float("inf"), 20, 10]),
            maxfev=5000,
        )
        T, B, logEC50, h = popt
        yh = _four_pl(log_x, *popt)
        ss_res = float(np.sum((y - yh) ** 2))
        ss_tot = float(np.sum((y - np.mean(y)) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        se = np.sqrt(np.diag(np.abs(pcov))) if np.all(np.isfinite(pcov)) else np.full(4, np.nan)
        ci_low = 10 ** (logEC50 - 1.96 * se[2]) if math.isfinite(se[2]) else float("nan")
        ci_high = 10 ** (logEC50 + 1.96 * se[2]) if math.isfinite(se[2]) else float("nan")
        return FourPLResult(
            Top=float(T), Bottom=float(B), logEC50=float(logEC50), HillSlope=float(h),
            EC50=10 ** logEC50, R2=r2, ssRes=ss_res, ssTot=ss_tot,
            nPoints=len(log_x), converged=True,
            seTop=float(se[0]), seBottom=float(se[1]),
            seLogEC50=float(se[2]), seHillSlope=float(se[3]),
            ci95Low=ci_low, ci95High=ci_high, method="scipy_curve_fit",
        )
    except Exception:
        return None


def _fit_pure_numpy(log_x, y):
    n = len(log_x)
    if n < 4:
        return None
    p = np.array(_initial_guess(log_x, y))
    lam = 1e-3

    def residuals(p):
        return y - _four_pl(log_x, *p)

    def jacobian(p):
        eps = 1e-7
        r0 = residuals(p)
        J = np.zeros((4, n))
        for j in range(4):
            pp = p.copy()
            pp[j] += eps
            re = residuals(pp)
            J[j] = (re - r0) / eps
        return J, r0

    converged = False
    iterations = 0
    for it in range(250):
        iterations = it + 1
        r = residuals(p)
        ss0 = float(np.sum(r * r))
        if ss0 < 1e-14:
            converged = True
            break
        J, _ = jacobian(p)
        JtJ = J @ J.T
        Jtr = J @ r
        JtJd = JtJ.copy()
        np.fill_diagonal(JtJd, np.diag(JtJd) * (1 + lam))
        try:
            delta = np.linalg.solve(JtJd, -Jtr)
        except np.linalg.LinAlgError:
            lam *= 5
            if lam > 1e10:
                break
            continue
        new_p = p + delta
        r_new = residuals(new_p)
        ss_new = float(np.sum(r_new * r_new))
        if ss_new < ss0:
            p = new_p
            lam = max(lam * 0.7, 1e-12)
            if it > 3 and np.max(np.abs(delta)) < 1e-10 and (ss0 - ss_new) / ss0 < 1e-10:
                converged = True
                break
        else:
            lam = min(lam * 5, 1e10)
            if lam > 1e10:
                break

    yh = _four_pl(log_x, *p)
    y_mean = float(np.mean(y))
    ss_res = float(np.sum((y - yh) ** 2))
    ss_tot = float(np.sum((y - y_mean) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    dof = max(1, n - 4)
    s2 = ss_res / dof
    J_final, _ = jacobian(p)
    JtJ_final = J_final @ J_final.T
    try:
        inv_JtJ = np.linalg.inv(JtJ_final)
        se = np.sqrt(np.abs(np.diag(inv_JtJ)) * s2)
    except np.linalg.LinAlgError:
        se = np.full(4, np.nan)
    T, B, logEC50, h = p
    se_logEC50 = se[2] if math.isfinite(se[2]) else float("nan")
    ci_low = 10 ** (logEC50 - 1.96 * se_logEC50) if math.isfinite(se_logEC50) else float("nan")
    ci_high = 10 ** (logEC50 + 1.96 * se_logEC50) if math.isfinite(se_logEC50) else float("nan")
    return FourPLResult(
        Top=float(T), Bottom=float(B), logEC50=float(logEC50), HillSlope=float(h),
        EC50=10 ** logEC50, R2=r2, ssRes=ss_res, ssTot=ss_tot,
        nPoints=n, converged=converged,
        seTop=float(se[0]), seBottom=float(se[1]),
        seLogEC50=se_logEC50, seHillSlope=float(se[3]),
        ci95Low=ci_low, ci95High=ci_high,
        iterations=iterations, method="pure_numpy_LM",
    )


def fit_4pl(log_conc, viability):
    """Fit a 4PL dose-response curve."""
    log_x = np.asarray(log_conc, dtype=float)
    y = np.asarray(viability, dtype=float)
    if HAS_SCIPY:
        result = _fit_scipy(log_x, y)
        if result is not None:
            return result
    return _fit_pure_numpy(log_x, y)


def predict_4pl(log_conc, top, bottom, logEC50, hill_slope):
    """Predict viability using fitted 4PL params."""
    return _four_pl(log_conc, top, bottom, logEC50, hill_slope)
