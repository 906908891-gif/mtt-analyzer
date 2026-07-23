"""Unit tests for statistical functions (mtt_analyzer.stats)."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import math
try:
    import pytest
except ImportError:
    pytest = None

from mtt_analyzer.stats import (
    mean, sd, sem, cv,
    welch_t_test, holm_bonferroni, significance_stars,
    incomplete_beta, _log_gamma,
)


def test_mean_basic():
    assert mean([1, 2, 3, 4, 5]) == 3.0
    assert mean([1.0]) == 1.0
    assert math.isnan(mean([]))


def test_sd_basic():
    assert sd([1, 2, 3, 4, 5]) == pytest.approx(1.5811388, rel=1e-4)
    assert sd([1.0]) == 0.0


def test_sem_basic():
    assert sem([1, 2, 3, 4, 5]) == pytest.approx(0.7071067, rel=1e-4)


def test_cv_basic():
    c = cv([1, 2, 3, 4, 5])
    assert c > 0
    assert math.isnan(cv([0, 0, 0]))


def test_welch_basic():
    a = [0.85, 0.87, 0.83, 0.86]
    b = [0.45, 0.42, 0.48, 0.43]
    result = welch_t_test(a, b)
    assert result["p"] < 0.001
    assert result["t"] > 0


def test_welch_same_groups():
    a = [0.85, 0.87, 0.83, 0.86]
    result = welch_t_test(a, a)
    assert result["p"] > 0.9
    assert abs(result["t"]) < 1e-9


def test_welch_degenerate():
    result = welch_t_test([1.0], [1.0])
    assert "p" in result


def test_holm_bonferroni_basic():
    p = [0.01, 0.04, 0.03, 0.005]
    adj = holm_bonferroni(p)
    assert all(a <= 1.0 for a in adj)
    assert all(a >= 0.0 for a in adj)


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


def test_log_gamma():
    assert _log_gamma(5) == pytest.approx(math.log(24), rel=1e-4)


def test_incomplete_beta_endpoints():
    assert incomplete_beta(1, 1, 0) == 0.0
    assert incomplete_beta(1, 1, 1) == 1.0


def test_incomplete_beta_uniform():
    assert incomplete_beta(1, 1, 0.5) == pytest.approx(0.5, rel=1e-6)
