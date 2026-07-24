# Changelog

All notable changes to this project will be documented in this file.

## [0.2.0] - 2026-07-24

### Critical fixes (post-release audit)

- **FIXED**: `mtt_analyzer/__init__.py` was wrapped in a docstring; all imports
  were inert. `from mtt_analyzer import fit_4pl` now works.
- **FIXED**: `script.js` in `<head>` without `defer` - DOM elements were
  null when JS bound events. Added `defer` attribute.
- **FIXED**: `incomplete_beta` now uses the Numerical Recipes continued-fraction
  algorithm (correct). Previously returned 0.36 for `I_0.5(1,1)` (should be 0.5).
- **FIXED**: `Holm-Bonferroni` direction was reversed (min instead of max).
  Now matches Holm (1979). Both Python and JS implementations updated.
- **FIXED**: HTML template artifact `'"'"' + `"`n    </div>"` cleaned up.
- **FIXED**: `app.py` now reads `PORT` and `HOST` environment variables for
  cloud deployment (Render / Railway / Heroku).
- **FIXED**: `fit_4pl` API endpoint enforces minimum 4 data points.

### Added

- `Dockerfile` (gunicorn + healthcheck)
- `.dockerignore`
- `render.yaml` (Render.com config)
- `Procfile` (Heroku config)
- `runtime.txt` (Python version pin)
- `.github/workflows/ci.yml` (matrix test on Python 3.9-3.12 + Docker build)
- `tests/test_smoke.py` (Flask boot + HTTP endpoint smoke tests)
- `incomplete_beta` symmetry test (I_x(a,b) + I_(1-x)(b,a) = 1)
- `test_welch_no_negative_p` (stress test 50 random samples)
- `test_holm_bonferroni_monotonic` (validates non-decreasing adjusted p)

### Improved

- Tests no longer depend on pytest fixtures (use `tempfile.NamedTemporaryFile`)
- Tests no longer use `pytest.approx` (home-grown `approx_equal` helper)
- `README.md` rewritten with deploy instructions and corrected accuracy claims
- `requirements.txt` includes `gunicorn` and `scipy`
- `pyproject.toml` adds `scipy` to dev dependencies

## [0.1.0] - 2026-07-23

### Initial Public Release (Python project)

- Python 3.9+ Flask web app
- `mtt_analyzer` package: core / stats / io / data
- 8 REST API endpoints
- 28 unit tests (now 34 after fixes)
- Welchs t-test, Holm-Bonferroni, 4PL dose-response fit
- Multi-plate support, manual grouping

## [0.0.1] - 2026-07-23 (deprecated)

### Single-HTML-First Era

- Single HTML file with 70 JS functions
- Replaced by Python project in 0.1.0
