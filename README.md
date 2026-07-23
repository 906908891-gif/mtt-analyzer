# MTT 实验数据分析器 (Python)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](./LICENSE)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)]()
[![Flask](https://img.shields.io/badge/Flask-2.0+-green.svg)]()
[![GitHub](https://img.shields.io/badge/GitHub-906908891--gif-blue.svg)](https://github.com/906908891-gif/mtt-analyzer)

一个用 **Python (Flask)** 写的 MTT 实验数据分析工具,提供网页界面 + REST API。

- 🧪 多板数据管理
- 📊 4PL 剂量响应拟合 + IC50 + 95% CI
- 📈 Welch's t-test + Holm-Bonferroni 多重比较校正
- 💾 localStorage 会话保存
- 🚀 一键部署到 Vercel / Render / Railway

## ✨ 核心特性

| 功能 | 描述 |
|------|------|
| 多板管理 | 同时处理多块 6/24/96 孔板,跨板合并数据 |
| 手动分组 | 可视化点击孔位 → 自定义实验组 / 对照组 |
| 排序 | 拖动 / 上下箭头重排,柱状图 X 轴自动跟随 |
| 统计检验 | Welch's t-test (不假设等方差) + Holm-Bonferroni |
| 4PL 拟合 | Levenberg-Marquardt 拟合,自动 IC50 + 95% CI |
| 图表 | 柱状图 / 热图 / 4PL 拟合曲线 |
| 导出 | CSV / PNG |
| API | REST endpoints 供外部脚本调用 |

## 🚀 快速开始

### 本地运行 (推荐)
```bash
git clone https://github.com/906908891-gif/mtt-analyzer.git
cd mtt-analyzer
pip install -r requirements.txt
python app.py
```
浏览器打开 http://localhost:5000。

### Docker
```dockerfile
FROM python:3.11-slim
WORKDIR /app
COPY . .
RUN pip install --no-cache-dir -r requirements.txt
EXPOSE 5000
CMD ["python", "app.py"]
```
```bash
docker build -t mtt-analyzer .
docker run -p 5000:5000 mtt-analyzer
```

### 编程使用
```python
from mtt_analyzer import fit_4pl

# 4PL 拟合
log_conc = [-6, -5, -4, -3, -2]  # log10 M
viability = [95.2, 78.3, 45.1, 12.4, 5.8]
result = fit_4pl(log_conc, viability)

print(f"EC50: {result.EC50:.2e} M")
print(f"95% CI: [{result.ci95Low:.2e}, {result.ci95High:.2e}]")
print(f"R2: {result.R2:.4f}")
print(f"Method: {result.method}")
```

```python
from mtt_analyzer.stats import welch_t_test, holm_bonferroni

# Welch t-test
t1 = [0.85, 0.87, 0.83]
t2 = [0.45, 0.42, 0.48]
result = welch_t_test(t1, t2)
print(f"p = {result['p']:.4f}")

# 多重校正
p_values = [0.001, 0.04, 0.12, 0.55]
adjusted = holm_bonferroni(p_values)
print(f"adjusted p = {adjusted}")
```

## 📂 项目结构

```
mtt-analyzer/
├── mtt_analyzer/          # 核心 Python 包
│   ├── __init__.py
│   ├── core.py            # 4PL 拟合 + IC50
│   ├── stats.py           # Welch t-test, Holm-Bonferroni
│   ├── io.py              # CSV / 浓度 / 孔位解析
│   └── data.py            # 数据模型 (Plate / Group / Session)
├── app.py                 # Flask 入口
├── templates/
│   └── index.html         # 主页面
├── static/
│   ├── style.css
│   └── script.js
├── tests/
│   ├── test_core.py       # 4PL 拟合单测
│   ├── test_stats.py      # 统计单测
│   └── test_io.py         # I/O 单测
├── examples/
│   └── sample.csv         # 示例数据
├── pyproject.toml         # Python 项目配置
├── requirements.txt
├── README.md
├── LICENSE
├── CHANGELOG.md
└── .gitignore
```

## 🧮 4PL 模型

```
y = Bottom + (Top - Bottom) / (1 + 10^((x - logEC50) * HillSlope))
```

最小二乘求解使用 **Levenberg-Marquardt** 算法 (`scipy.optimize.curve_fit`),无 scipy 时自动降级到纯 numpy 实现。

## 🔌 REST API

| 端点 | 方法 | 功能 |
|------|------|------|
| `/api/health` | GET | 健康检查 |
| `/api/parse-wells` | POST | 解析 `A1-A6, B1-B3` 格式 |
| `/api/parse-concentration` | POST | 解析 `1uM` → `1e-6` |
| `/api/welch-ttest` | POST | Welch's t-test |
| `/api/holm-bonferroni` | POST | 多重校正 |
| `/api/fit-4pl` | POST | 4PL 拟合 + IC50 + 95% CI |
| `/api/predict-4pl` | POST | 用拟合参数预测 |
| `/api/compute-stats` | POST | 批量统计 |

## 🧪 运行测试

```bash
pip install -e ".[dev]"
pytest tests/
```

## 📝 License

[MIT](./LICENSE) (c) 2026 906908891-gif
