# scripts/figures —— 十张图的唯一绘图实现

`manuscript_figures.py` 是**作者端原始绘图代码的移植版**：把各轮（15/34/48/49/51/52，以及
legacy 435 R 脚本）里的绘图代码段逐段搬进本仓库，只做两件事——

1. 输入路径改为相对路径（`data/derived`、`data/upstream`）；
2. 数值列**按列名**取用（不使用位置索引 `iloc[:, n]`）。

绘图参数（面板划分、配色、标签、统计注脚、坐标范围）与原代码保持一致，因此产出图件与稿件
`manuscript_v6/assets/` 一一对应。

## 逐图来源索引

| 函数 | 原始代码 | 备注 |
|---|---|---|
| `fig1_cohort_strategies` | 51 `build_figures.py` Figure 1 段 | a 队列计数 / b 门组成 / c 模块 breadth 热图（含数值标注） |
| `fig2_lineage_mosaic` | 51 `restore_tree_figure.py` 整文件 | 原始细菌/古菌树 + 谱系标签 + breadth 热图 |
| `fig3_portfolio_uncertainty` | 48 `plot_diagnostics.py` 整文件 | A 配对不确定区间 / B 模块分解 / C 面板敏感性 + 原始注脚 |
| `fig4_molecular_configurations` | 51 `build_figures.py` Figure 4 段 | a 共检率（系统 × 来源，含分母）/ b family 层 partial R² / c 邻近代理 |
| `fig5_dryland_context` | 51 `build_figures.py` Figure 5 段 | a 标记集合 / b 10 kb 共定位 / c 谱系支持 / d 基因箭头位点 |
| `figs1_framework` | 34 `build_figures.py`（framework） | 概念图，无数据输入 |
| `figs2_bayesian_profiles` | 34 `build_figures.py`（bayesian） | 5 行后验中位数/区间 |
| `figs3_null_details` | 51 `build_figures.py` KO-null 段 | 仅 within-source 零模型；c 面板含 global 对比 |
| `figs4_independent_occurrence` | legacy `435.redraw_corrected_WGS_figure_20260927.R` | 按 (A\|B)/C 版式等价重绘 |
| `fig6_dpfq008` | 52 `restore_candidate_figure.py` | 原始矢量 + 已批准标签修正链 |

每个函数签名统一为 `plot(root: Path, out: Path) -> list[str]`，由 `bin/nee_recon.py plot <label>` 调用。

## 不做什么

- 不复制稿件里的成品 PNG 充当"复现"（Fig6 是算法化的原始矢量 + 标签修正链，属特例且已说明）。
- 不重跑任何重型分析、不改动任何科学阈值。
- 不引入新依赖：仅用 `numpy / pandas / matplotlib`（Fig6 的栅格化另需 `rsvg-convert` 或 `cairosvg`）。
