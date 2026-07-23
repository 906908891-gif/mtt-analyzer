"""MTT Analyzer - Flask web app entry point.

Run locally:
    python app.py

Then open http://localhost:5000 in a browser.
"""

from __future__ import annotations

import argparse

from flask import Flask, render_template, request, jsonify

from mtt_analyzer import fit_4pl, predict_4pl
from mtt_analyzer.stats import welch_t_test, holm_bonferroni, mean, sd, sem, cv, significance_stars
from mtt_analyzer.io import parse_concentration, expand_wells


def create_app():
    app = Flask(__name__, static_folder="static", template_folder="templates")

    @app.route("/")
    def index():
        return render_template("index.html")

    @app.route("/api/health")
    def health():
        return jsonify({"status": "ok"})

    @app.route("/api/parse-wells", methods=["POST"])
    def api_parse_wells():
        data = request.get_json(silent=True) or {}
        pattern = data.get("pattern", "")
        rows = int(data.get("rows", 8))
        cols = int(data.get("cols", 12))
        try:
            wells = expand_wells(pattern, (rows, cols))
        except Exception as e:
            return jsonify({"error": str(e)}), 400
        return jsonify({"wells": wells, "count": len(wells)})

    @app.route("/api/parse-concentration", methods=["POST"])
    def api_parse_concentration():
        data = request.get_json(silent=True) or {}
        text = data.get("text", "")
        val = parse_concentration(text)
        return jsonify({"value": val})

    @app.route("/api/welch-ttest", methods=["POST"])
    def api_welch_ttest():
        data = request.get_json(silent=True) or {}
        a = data.get("a", [])
        b = data.get("b", [])
        result = welch_t_test(a, b)
        return jsonify(result)

    @app.route("/api/holm-bonferroni", methods=["POST"])
    def api_holm_bonferroni():
        data = request.get_json(silent=True) or {}
        p_values = data.get("p_values", [])
        adjusted = holm_bonferroni(p_values)
        return jsonify({"adjusted": adjusted})

    @app.route("/api/fit-4pl", methods=["POST"])
    def api_fit_4pl():
        data = request.get_json(silent=True) or {}
        log_conc = data.get("log_conc", [])
        viability = data.get("viability", [])
        if not log_conc or not viability or len(log_conc) != len(viability):
            return jsonify({"error": "log_conc and viability must be non-empty and equal length"}), 400
        result = fit_4pl(log_conc, viability)
        if result is None:
            return jsonify({"error": "fit failed (need at least 4 points)"}), 400
        return jsonify(result.to_dict())

    @app.route("/api/predict-4pl", methods=["POST"])
    def api_predict_4pl():
        data = request.get_json(silent=True) or {}
        log_conc = data.get("log_conc", [])
        params = data.get("params", {})
        if not all(k in params for k in ("top", "bottom", "logEC50", "hill_slope")):
            return jsonify({"error": "missing params"}), 400
        y = predict_4pl(
            log_conc,
            float(params["top"]), float(params["bottom"]),
            float(params["logEC50"]), float(params["hill_slope"]),
        )
        return jsonify({"predicted": [float(v) for v in y]})

    @app.route("/api/compute-stats", methods=["POST"])
    def api_compute_stats():
        data = request.get_json(silent=True) or {}
        groups = data.get("groups", [])
        control_id = data.get("control_id")
        if not isinstance(groups, list):
            return jsonify({"error": "groups must be a list"}), 400
        control = next((g for g in groups if g.get("id") == control_id), None)
        if not control:
            return jsonify({"error": "control group not found"}), 400
        control_values = [float(v) for v in control.get("values", [])]
        if len(control_values) < 2:
            return jsonify({"error": "control group needs >= 2 values"}), 400

        results = []
        for g in groups:
            values = [float(v) for v in g.get("values", [])]
            if not values:
                continue
            results.append({
                "id": g.get("id"),
                "name": g.get("name", ""),
                "color": g.get("color", "#0891b2"),
                "n": len(values),
                "mean": mean(values),
                "sd": sd(values),
                "sem": sem(values),
                "cv": cv(values),
                "values": values,
            })

        ctrl_stats = next((r for r in results if r["id"] == control_id), None)
        if not ctrl_stats:
            return jsonify({"error": "control stats missing"}), 400
        comp_stats = [r for r in results if r["id"] != control_id]
        p_raw = [welch_t_test(r["values"], control_values)["p"] for r in comp_stats]
        p_adj = holm_bonferroni(p_raw)
        for r, p, pa in zip(comp_stats, p_raw, p_adj):
            r["pValue"] = p
            r["pAdjusted"] = pa
            r["significance"] = significance_stars(pa)
            r["viability"] = (r["mean"] / ctrl_stats["mean"]) * 100.0
        ctrl_stats["viability"] = 100.0
        ctrl_stats["pValue"] = None
        ctrl_stats["pAdjusted"] = None
        ctrl_stats["significance"] = ""
        return jsonify({"groups": results, "control": ctrl_stats})

    return app


def main():
    parser = argparse.ArgumentParser(description="MTT Analyzer Flask web app")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    app = create_app()
    print(f"  MTT Analyzer: http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
