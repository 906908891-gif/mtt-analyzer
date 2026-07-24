"""Smoke tests for Flask app and its HTTP endpoints."""

import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest

from app import create_app


@pytest.fixture
def client():
    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as c:
        yield c


def test_health(client):
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.get_json()["status"] == "ok"


def test_index_renders(client):
    r = client.get("/")
    assert r.status_code == 200
    body = r.get_data(as_text=True)
    assert "<!DOCTYPE html>" in body
    assert "MTT" in body
    assert "url_for" in body  # Flask template rendered


def test_index_has_defer_script(client):
    """script.js must be loaded with defer to avoid DOM-not-ready bugs."""
    r = client.get("/")
    body = r.get_data(as_text=True)
    assert chr(60) + "script src=" + chr(34) in body
    assert "defer" in body, "script.js should have defer attribute"


def test_parse_wells_endpoint(client):
    r = client.post("/api/parse-wells", json={"pattern": "A1-A3, B1"})
    assert r.status_code == 200
    data = r.get_json()
    assert "wells" in data
    assert data["count"] == 4
    assert "A3" in data["wells"]


def test_parse_concentration_endpoint(client):
    r = client.post("/api/parse-concentration", json={"text": "10uM"})
    assert r.status_code == 200
    data = r.get_json()
    assert data["value"] == 1e-5


def test_welch_ttest_endpoint(client):
    r = client.post("/api/welch-ttest", json={
        "a": [0.85, 0.87, 0.83],
        "b": [0.45, 0.42, 0.48],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "p" in data
    assert data["p"] < 0.001
    assert data["p"] >= 0.0


def test_holm_bonferroni_endpoint(client):
    r = client.post("/api/holm-bonferroni", json={"p_values": [0.01, 0.04, 0.03, 0.005]})
    assert r.status_code == 200
    data = r.get_json()
    adj = data["adjusted"]
    # Expected: [0.03, 0.06, 0.06, 0.02]
    expected = [0.03, 0.06, 0.06, 0.02]
    for got, exp in zip(adj, expected):
        assert abs(got - exp) < 1e-4, (got, exp)


def test_fit_4pl_endpoint(client):
    r = client.post("/api/fit-4pl", json={
        "log_conc": [-5, -4, -3, -2, -1],
        "viability": [85, 60, 35, 15, 5],
    })
    assert r.status_code == 200
    data = r.get_json()
    assert "EC50" in data
    assert "R2" in data
    assert "Top" in data
    assert data["R2"] > 0.95


def test_fit_4pl_too_few(client):
    r = client.post("/api/fit-4pl", json={
        "log_conc": [-5, -4, -3],
        "viability": [80, 50, 20],
    })
    assert r.status_code == 400


def test_compute_stats_endpoint(client):
    r = client.post("/api/compute-stats", json={
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
    # g2 has p-value
    g2 = next(g for g in data["groups"] if g["id"] == "g2")
    assert g2["pValue"] < 0.001
    assert g2["pAdjusted"] >= g2["pValue"]  # Holm only increases p
