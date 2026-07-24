"""Unit tests for statistical functions (mtt_analyzer.stats)."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import math
import numpy as np

from mtt_analyzer.stats import (
    mean, sd, sem, cv,
    welch_t_test, holm_bonferroni, significance_stars,
    incomplete_beta,
)


def approx_equal(a, b, rel=1e-6, abs_tol=1e-9):
    """Compare two floats with relative and absolute tolerance."""
    if abs(a - b) < abs_tol:
        return True
    return abs(a - b) / max(abs(a), abs(b), 1e-30) < rel


def test_mean_basic():
    assert mean([1, 2, 3, 4, 5]) == 3.0
    assert mean([1.0]) == 1.0
    assert math.isnan(mean([]))


def test_sd_basic():
    assert approx_equal(sd([1, 2, 3, 4, 5]), 1.5811388)
    assert sd([1.0]) == 0.0


def test_sem_basic():
    assert approx_equal(sem([1, 2, 3, 4, 5]), 0.7071067)


def test_cv_basic():
    c = cv([1, 2, 3, 4, 5])
    assert c > 0
    assert math.isnan(cv([0, 0, 0]))


def test_welch_basic():
    a = [0.85, 0.87, 0.83, 0.86]
    b = [0.45, 0.42, 0.48, 0.43]
    result = welch_t_test(a, b)
    assert result["p"] < 0.001
    assert 0.0 <= result["p"] <= 1.0
    assert result["t"] > 0


def test_welch_same_groups():
    a = [0.85, 0.87, 0.83, 0.86]
    result = welch_t_test(a, a)
    assert result["p"] > 0.9
    assert abs(result["t"]) < 1e-9


def test_welch_degenerate():
    result = welch_t_test([1.0], [1.0])
    assert "p" in result
    assert result["p"] == 1.0


def test_welch_no_negative_p():
    """Stress test: p must be in [0, 1] for many random samples."""
    import random
    random.seed(42)
    for _ in range(50):
        a = [random.gauss(0.5, 0.1) for _ in range(5)]
        b = [random.gauss(0.3, 0.1) for _ in range(5)]
        r = welch_t_test(a, b)
        assert 0.0 <= r["p"] <= 1.0, "p out of range"


def test_holm_bonferroni_basic():
    """Holm-Bonferroni reference values (p=[0.01, 0.04, 0.03, 0.005])."""
    p = [0.01, 0.04, 0.03, 0.005]
    adj = holm_bonferroni(p)
    expected = [0.03, 0.06, 0.06, 0.02]
    for got, exp in zip(adj, expected):
        assert approx_equal(got, exp, rel=1e-4), (got, exp)


def test_holm_bonferroni_monotonic():
    """Adjusted p-values (sorted by original) must be non-decreasing."""
    p = [0.001, 0.5, 0.3, 0.05, 0.8]
    adj = holm_bonferroni(p)
    sorted_adj = [adj[i] for i in np.argsort(p)]
    for i in range(len(sorted_adj) - 1):
        assert sorted_adj[i] <= sorted_adj[i + 1] + 1e-9, sorted_adj


def test_holm_bonferroni_single():
    assert holm_bonferroni([0.05]) == [0.05]


def test_holm_bonferroni_empty():
    assert holm_bonferroni([]) == []


def test_significance_stars():
    assert significance_stars(0.0001) == "***"
    assert significance_stars(0.005) == "**"
    assert significance_stars(0.02) == "*"
    assert significance_stars(0.07) == "."
    assert significance_stars(0.5) == "ns"
    assert significance_stars(None) == ""
    assert significance_stars(float("nan")) == ""


def test_incomplete_beta_endpoints():
    assert incomplete_beta(1, 1, 0) == 0.0
    assert incomplete_beta(1, 1, 1) == 1.0
    assert incomplete_beta(2, 3, 0) == 0.0
    assert incomplete_beta(2, 3, 1) == 1.0


def test_incomplete_beta_uniform():
    """I_0.5(1, 1) = 0.5 (uniform on [0, 1])."""
    assert approx_equal(incomplete_beta(1, 1, 0.5), 0.5)


def test_incomplete_beta_reference():
    """Known reference values computed from scipy.special.betainc."""
    # I_0.4(2, 3) = 0.5248
    assert approx_equal(incomplete_beta(2, 3, 0.4), 0.5248, rel=1e-3)
    # I_0.5(0.5, 0.5) = 0.5 (symmetric)
    assert approx_equal(incomplete_beta(0.5, 0.5, 0.5), 0.5)
    # I_0.7(3, 2) = 0.652 (computed via numerical integration)
    assert approx_equal(incomplete_beta(3, 2, 0.7), 0.652, rel=1e-2)


def test_incomplete_beta_symmetry():
    """I_x(a, b) + I_(1-x)(b, a) = 1 (symmetry)."""
    for a, b in [(1, 1), (2, 3), (3, 2), (0.5, 0.5)]:
        for x in [0.1, 0.3, 0.5, 0.7, 0.9]:
            s = incomplete_beta(a, b, x) + incomplete_beta(b, a, 1 - x)
            assert approx_equal(s, 1.0, rel=1e-4), (a, b, x, s)
