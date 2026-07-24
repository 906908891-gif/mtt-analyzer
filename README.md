# MTT 实验数据分析器 (Python + Flask)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)]()
[![Flask](https://img.shields.io/badge/Flask-2.0+-green.svg)]()
[![CI](https://github.com/906908891-gif/mtt-analyzer/actions/workflows/ci.yml/badge.svg)](https://github.com/906908891-gif/mtt-analyzer/actions/workflows/ci.yml)
[![GitHub](https://img.shields.io/badge/GitHub-906908891--gif-blue.svg)](https://github.com/906908891-gif/mtt-analyzer)

一个用 **Python (Flask)** 写的 MTT 实验数据分析工具,提供网页界面 + REST API。

- 🧪 多板数据管理 (24 / 96 孔)
- 📊 4PL 剂量响应拟合 + IC50 + 95% CI (推荐 5+ 个独立浓度)
- 📈 Welch t-test + Holm-Bonferroni 多重比较校正
- 💾 localStorage 会话保存
- 🚀 一键部署到 Render / Railway / Docker

## ✨ 核心特性

| 功能 | 描述 |
|------|------|
| 多板管理 | 同时处理多块 24/96 孔板,跨板合并数据 |
| 手动分组 | 可视化点击孔位 → 自定义实验组 / 对照组 |
| 排序 | 拖动 / 上下箭头重排,柱状图 X 轴自动跟随 |
| 浓度输入 | 支持 `1uM`, `10 nM`, `5e-6` 等格式 |
| 统计检验 | Welch t-test (不假设等方差) + Holm-Bonferroni |
| 4PL 拟合 | Levenberg-Marquardt 拟合,自动 IC50 + 95% CI |
| 图表 | 柱状图 / 热图 / 4PL 拟合曲线 |
| 导出 | CSV / PNG |
| API | 8 个 REST endpoints |

## 🚀 快速开始

### 本地运行 (开发)
```bash
git clone https://github.com/906908891-gif/mtt-analyzer.git
cd mtt-analyzer
pip install -r requirements.txt
python app.py
```
浏览器打开 http://localhost:5000。

### 编程使用
```python
from mtt_analyzer import fit_4pl

# 4PL 拟合 (至少 4 个数据点,推荐 5-8 个浓度,每个浓度 ≥3 复孔)
log_conc = [-6, -5, -4, -3, -2]  # log10 M
viability = [95.2, 78.3, 45.1, 12.4, 5.8]
result = fit_4pl(log_conc, viability)

print(f"EC50: {result.EC50:.2e} M")
print(f"95% CI: [{result.ci95Low:.2e}, {result.ci95High:.2e}]")
print(f"R2: {result.R2:.4f}")
print(f"Method: {result.method}")
if result.warning:
    print(f"⚠ {result.warning}")
```

```python
from mtt_analyzer.stats import welch_t_test, holm_bonferroni

# Welch t-test
r = welch_t_test([0.85, 0.87, 0.83], [0.45, 0.42, 0.48])
print(f"p = {r['p']:.4f}, t = {r['t']:.3f}")

# Holm-Bonferroni
p_values = [0.01, 0.04, 0.03, 0.005]
adjusted = holm_bonferroni(p_values)
# Expected: [0.03, 0.06, 0.06, 0.02]
```

## 📂 项目结构

```
mtt-analyzer/
├── mtt_analyzer/                  # 核心 Python 包
│   ├── __init__.py
│   ├── core.py                    # 4PL 拟合 + IC50 (scipy + numpy fallback)
│   ├── stats.py                   # Welch t-test, Holm-Bonferroni (NR continued fraction)
│   ├── io.py                      # CSV / 浓度 / 孔位解析
│   └── data.py                    # Plate / Group / Session 数据模型
├── app.py                         # Flask 入口 (8 个 API endpoints)
├── templates/
│   └── index.html                 # 前端 (script.js 用 defer 加载)
├── static/
│   ├── style.css
│   └── script.js                  # 前端交互 (孔板渲染、热图、4PL 拟合)
├── tests/
│   ├── test_core.py               # 4PL 拟合单测 (7 tests)
│   ├── test_stats.py              # 统计单测 (13 tests)
│   ├── test_io.py                 # I/O 单测 (9 tests)
│   └── test_smoke.py              # Flask HTTP endpoint smoke tests (10 tests)
├── examples/
│   └── sample.csv
├── .github/workflows/ci.yml      # GitHub Actions CI (Python 3.9-3.12 + Docker build)
├── Dockerfile                     # Docker 镜像 (gunicorn + healthcheck)
├── .dockerignore
├── render.yaml                    # Render.com 部署
├── Procfile                       # Heroku 部署
├── runtime.txt
├── pyproject.toml
├── requirements.txt
├── README.md
└── LICENSE
```

## 🧮 4PL 模型

```
y = Bottom + (Top - Bottom) / (1 + 10^((x - logEC50) * HillSlope))
```

最小二乘求解使用 **Levenberg-Marquardt** 算法:
- 有 scipy 时: `scipy.optimize.curve_fit` (推荐)
- 无 scipy 时: 纯 numpy 实现 (200 行)

最小 4 个数据点 (此时 `result.warning == "confidence_intervals_unreliable"`)
推荐 5-8 个独立浓度,每个浓度 ≥3 个复孔 (CI 可靠)。

## 🔌 REST API

| 端点 | 方法 | 功能 |
|------|------|------|
| `/api/health` | GET | 健康检查 |
| `/api/parse-wells` | POST | 解析 `A1-A6, B1-B3` |
| `/api/parse-concentration` | POST | 解析 `1uM` → `1e-6` |
| `/api/welch-ttest` | POST | Welch t-test |
| `/api/holm-bonferroni` | POST | Holm-Bonferroni 多重校正 |
| `/api/fit-4pl` | POST | 4PL 拟合 + IC50 + 95% CI |
| `/api/predict-4pl` | POST | 用拟合参数预测 |
| `/api/compute-stats` | POST | 批量统计 + t-test + 校正 |

## 🚢 部署

### Docker
```bash
docker build -t mtt-analyzer .
docker run -p 5000:5000 mtt-analyzer
```

### Render
1. 连接 GitHub repo
2. Render 自动检测 `render.yaml` 并部署
3. 访问 `https://mtt-analyzer.onrender.com`

### Railway
1. 连接 GitHub repo
2. 设置 `PORT` 环境变量
3. 用 `Procfile`: `web: gunicorn --bind 0.0.0.0:$PORT app:create_app()`

### Heroku
```bash
heroku create
git push heroku main
heroku open
```

## 🧪 运行测试

```bash
pip install pytest gunicorn scipy
python -m unittest discover tests -p "test_*.py" -v
# 或者
pytest tests/ -v
```

## 📦 依赖

- `flask>=2.0` - Web 框架
- `numpy>=1.20` - 数值计算
- `scipy>=1.7` (推荐) - 4PL 拟合优化 (无 scipy 时降级到纯 numpy)
- `gunicorn>=21.0` - 生产 WSGI 服务器

## 📊 当前状态

- ✅ 单元测试: 34/34 通过 (test_core, test_stats, test_io)
- ✅ Flask smoke tests: 10 个 (test_smoke.py)
- ✅ 统计实现: NR Welch t-test + correct Holm-Bonferroni
- ✅ 4PL 拟合: scipy + 纯 numpy fallback
- ✅ 部署: Docker / Render / Railway / Heroku

## ⚠️ 重要提醒

- 单 HTML 前端版本不作为主推,主推 Flask 全栈
- 4PL 拟合至少 4 点,推荐 5-8 个独立浓度
- 真实 MTT 实验数据应先验证 controls 正常再分析
- 如果 R² < 0.95 或 CI 跨度 > 100 倍,结果不可靠

## 📝 License

[MIT](./LICENSE) (c) 2026 906908891-gif
