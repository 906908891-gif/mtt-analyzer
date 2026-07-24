"""Smoke tests for Flask app and its HTTP endpoints."""

import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from app import create_app


def _client():
    app = create_app()
    app.config["TESTING"] = True
    return app.test_client()


def test_health():
    c = _client()
    r = c.get("/api/health")
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"


def test_index_renders():
    c = _client()
    r = c.get("/")
    assert r.status_code == 200
    body = r.get_data(as_text=True)
    assert "<!DOCTYPE html>" in body
    assert "MTT" in body
    # Jinja should have rendered url_for -> /static/...
    assert "/static/style.css" in body
    assert "/static/script.js" in body
    # No unrendered Jinja placeholders
    assert "{{" not in body
    assert "{%" not in body


def test_index_has_defer_script():
    """script.js must be loaded with defer to avoid DOM-not-ready bugs."""
    c = _client()
    r = c.get("/")
    body = r.get_data(as_text=True)
    assert chr(60) + "script src=" + chr(34) in body
    assert "defer" in body


def test_parse_wells_endpoint():
    c = _client()
    r = c.post("/api/parse-wells", json={"pattern": "A1-A3, B1"})
    assert r.status_code == 200
    data = r.get_json()
    assert "wells" in data
    assert data["count"] == 4
    assert "A3" in data["wells"]


def test_parse_concentration_endpoint():
    c = _client()
    r = c.post("/api/parse-concentration", json={"text": "10uM"})
    assert r.status_code == 200
    data = r.get_json()
    # Floating point: just check it's close to 1e-5
    assert math.isclose(data["value"], 1e-5, rel_tol=1e-12, abs_tol=1e-15)


def test_welch_ttest_endpoint():
    c = _client()
    r = c.post("/api/welch-ttest", json={
        "a": [0.85, 0.87, 0.83],
        "b": [0.45, 0.42, 0.48],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "p" in data
    assert data["p"] < 0.001
    assert data["p"] >= 0.0


def test_holm_bonferroni_endpoint():
    c = _client()
    r = c.post("/api/holm-bonferroni", json={"p_values": [0.01, 0.04, 0.03, 0.005]})
    assert r.status_code == 200
    data = r.get_json()
    adj = data["adjusted"]
    expected = [0.03, 0.06, 0.06, 0.02]
    for got, exp in zip(adj, expected):
        assert abs(got - exp) < 1e-4


def test_fit_4pl_endpoint():
    c = _client()
    r = c.post("/api/fit-4pl", json={
        "log_conc": [-5, -4, -3, -2, -1],
        "viability": [85, 60, 35, 15, 5],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "EC50" in data
    assert "R2" in data
    assert "Top" in data
    assert data["R2"] > 0.95


def test_fit_4pl_too_few():
    c = _client()
    r = c.post("/api/fit-4pl", json={
        "log_conc": [-5, -4, -3],
        "viability": [80, 50, 20],
    })
    assert r.status_code == 400


def test_fit_4pl_four_points_warns():
    """With 4 points fit works but CI may be unreliable."""
    c = _client()
    r = c.post("/api/fit-4pl", json={
        "log_conc": [-5, -4, -3, -2],
        "viability": [80, 50, 20, 5],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "EC50" in data
    assert "warning" in data
    assert data["warning"] == "confidence_intervals_unreliable"


def test_compute_stats_endpoint():
    c = _client()
    r = c.post("/api/compute-stats", json={
        "control_id": "g1",
        "groups": [
            {"id": "g1", "name": "Control", "values": [0.85, 0.87, 0.83]},
            {"id": "g2", "name": "Trt", "values": [0.45, 0.42, 0.48]},
        ],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "groups" in data
    assert len(data["groups"]) == 2
    g2 = next(g for g in data["groups"] if g["id"] == "g2")
    assert g2["pValue"] < 0.001
    assert g2["pAdjusted"] >= g2["pValue"]


def test_compute_stats_handles_single_replicate():
    c = _client()
    r = c.post("/api/compute-stats", json={
        "control_id": "g1",
        "groups": [
            {"id": "g1", "name": "Control", "values": [0.85, 0.87, 0.83, 0.86]},
            {"id": "g2", "name": "TrtA", "values": [0.45, 0.42, 0.48]},
            {"id": "g3", "name": "Single", "values": [0.5]},
            {"id": "g4", "name": "TrtC", "values": [0.30, 0.32, 0.28]},
        ],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert len(data["groups"]) == 4
    # g3 single-replicate -> marked 样本不足
    g3 = next(g for g in data["groups"] if g["id"] == "g3")
    assert g3["note"] == "样本不足"
    assert g3["pValue"] is None
    assert g3["pAdjusted"] is None
    assert g3["significance"] == "N/A"
    # g2 and g4 (n>=2) get p-values
    g2 = next(g for g in data["groups"] if g["id"] == "g2")
    assert g2["pValue"] is not None
    assert g2["pAdjusted"] is not None
    g4 = next(g for g in data["groups"] if g["id"] == "g4")
    assert g4["pValue"] is not None
    # Holm-adjusted p should only consider 2 groups (g2, g4), not g3
    assert g2["pAdjusted"] >= g2["pValue"]


def test_compute_stats_control_too_small():
    c = _client()
    r = c.post("/api/compute-stats", json={
        "control_id": "g1",
        "groups": [
            {"id": "g1", "name": "Control", "values": [0.85]},
            {"id": "g2", "name": "Trt", "values": [0.45, 0.42, 0.48]},
        ],
    })
    assert r.status_code == 400
