"""Unit tests for 4PL fit (mtt_analyzer.core)."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import math
try:
    import pytest
except ImportError:
    pytest = None
import numpy as np

from mtt_analyzer.core import fit_4pl, predict_4pl, FourPLResult


def _make_4pl_data(log_conc, top=100.0, bottom=0.0, logEC50=-4.0, hill=1.0, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(log_conc, dtype=float)
    y = bottom + (top - bottom) / (1.0 + np.power(10.0, (x - logEC50) * hill))
    if noise > 0:
        y += rng.normal(0, noise, size=y.shape)
    return x, y


def test_fit_perfect_data():
    log_conc = [-7, -6, -5, -4, -3, -2]
    _, y = _make_4pl_data(log_conc, top=100.0, bottom=0.0, logEC50=-4.0, hill=1.0)
    result = fit_4pl(log_conc, y)
    assert result is not None
    assert result.converged
    assert abs(result.EC50 - 1e-4) < 1e-6
    assert abs(result.Top - 100.0) < 0.01
    assert abs(result.Bottom - 0.0) < 0.01
    assert abs(result.HillSlope - 1.0) < 0.01
    assert result.R2 > 0.999


def test_fit_too_few_points():
    result = fit_4pl([-5, -4, -3], [50, 40, 30])
    assert result is None


def test_fit_with_noise():
    log_conc = [-7, -6, -5, -4, -3, -2, -1]
    _, y = _make_4pl_data(log_conc, top=100.0, bottom=5.0, logEC50=-4.5, hill=1.2, noise=2.0, seed=42)
    result = fit_4pl(log_conc, y)
    assert result is not None
    assert abs(result.EC50 - 10 ** -4.5) / 10 ** -4.5 < 0.1
    assert result.R2 > 0.95


def test_predict_4pl():
    log_conc = [-5, -4, -3]
    y = np.array([80.0, 50.0, 20.0])
    result = fit_4pl(log_conc, y)
    yh = predict_4pl(log_conc, result.Top, result.Bottom, result.logEC50, result.HillSlope)
    assert np.allclose(yh, y, atol=1e-3)


def test_result_to_dict():
    log_conc = [-5, -4, -3, -2]
    y = [80.0, 50.0, 20.0, 5.0]
    result = fit_4pl(log_conc, y)
    d = result.to_dict()
    assert "EC50" in d
    assert "R2" in d
    assert "method" in d
    assert isinstance(d["EC50"], float)


def test_four_pl_result_fields():
    log_conc = [-5, -4, -3, -2]
    y = [80.0, 50.0, 20.0, 5.0]
    result = fit_4pl(log_conc, y)
    assert isinstance(result, FourPLResult)
    assert math.isfinite(result.EC50)
    assert 0.0 <= result.R2 <= 1.0
