"""Unit tests for I/O helpers (mtt_analyzer.io)."""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import os
try:
    import pytest
except ImportError:
    pytest = None
import tempfile

from mtt_analyzer.io import parse_concentration, parse_csv, expand_wells


def test_parse_concentration_pure_number():
    assert parse_concentration("1e-6") == 1e-6
    assert parse_concentration("0.5") == 0.5
    assert parse_concentration("100") == 100.0


def test_parse_concentration_si():
    assert parse_concentration("1uM") == 1e-6
    assert parse_concentration("10 nM") == 1e-8
    assert parse_concentration("1mM") == 1e-3
    assert parse_concentration("5M") == 5.0
    assert parse_concentration("1pM") == 1e-12


def test_parse_concentration_invalid():
    assert parse_concentration("") is None
    assert parse_concentration(None) is None
    assert parse_concentration("abc") is None


def test_parse_csv(tmp_path):
    f = tmp_path / "data.csv"
    csv_content = "plate,well,group\nA1,Control\nB1,Treatment\n"
    f.write_text(csv_content, encoding="utf-8")
    rows = parse_csv(str(f))
    assert len(rows) == 2
    assert rows[0]["plate"] == "A1"
    assert rows[0]["group"] == "Control"
    assert rows[1]["well"] == "B1"


def test_expand_wells_single():
    assert expand_wells("A1") == ["A1"]
    assert expand_wells("A1, A2") == ["A1", "A2"]


def test_expand_wells_range():
    assert expand_wells("A1-A3") == ["A1", "A2", "A3"]
    assert expand_wells("A1-A3, B1") == ["A1", "A2", "A3", "B1"]


def test_expand_wells_bracket():
    result = expand_wells("[A-C]5")
    assert "A5" in result
    assert "B5" in result
    assert "C5" in result
    assert len(result) == 3


def test_expand_wells_dedup():
    result = expand_wells("A1, A1, A2")
    assert result == ["A1", "A2"]


def test_expand_wells_empty():
    assert expand_wells("") == []
    assert expand_wells("   ") == []
