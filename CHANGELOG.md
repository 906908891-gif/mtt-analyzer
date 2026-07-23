# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased] - 2026-07-23

### Changed
- **Project type**: Migrated from single-HTML-file to Python Flask project
- **Analysis backend**: Ported 4PL fit, Welch t-test, Holm-Bonferroni from JavaScript to Python
- **Project structure**: Now a proper Python package with `pyproject.toml`

### Added
- `mtt_analyzer/` Python package with sub-modules:
  - `core.py` - 4PL fit (scipy + pure-numpy fallback)
  - `stats.py` - Welch t-test, Holm-Bonferroni, basic stats
  - `io.py` - CSV / concentration / well-pattern parsing
  - `data.py` - Data classes for Plate / Group / Session
- `app.py` - Flask web app with REST API endpoints
- `tests/` - Unit tests for core, stats, io
- `pyproject.toml` and `requirements.txt` - Python packaging
- `examples/sample.csv` - Sample data for quick testing
- REST API endpoints: `/api/parse-wells`, `/api/fit-4pl`, `/api/welch-ttest`, etc.

## [0.0.1] - 2026-07-23 (deprecated)

### Single-HTML-First Era
- 单 HTML 文件实现, 70 个 JS 函数
- 4PL 拟合用 Levenberg-Marquardt (JS 版)
- Welch t-test + Holm-Bonferroni in JS
- GitHub Pages 部署
