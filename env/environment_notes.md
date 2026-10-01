# 计算环境记录（env/environment_notes.md）

本文件记录复现本论文 v6 图件与核心表所需的**真实**软件环境与来源。所有版本号均来自项目内既有的
`requirements.lock.txt`（15 轮 Bayesian 分析，真实 pip 冻结）、52 轮 `manuscript_v6.md` 的 Methods 段，
以及本仓库脚本源码中的引用；未联网核对者一律标注 `(未联网核对)`。

## Python 环境（图件与 WGS 下游）

`env/requirements.lock.txt` 由 `15_bayesian_environment_strategies_20260927/requirements.lock.txt` 原样复制
（真实 `pip freeze` 输出）。关键包：

| 包 | 版本 |
|---|---|
| python | 3.12（venv 名为 `.venv-nee-bayes`，位于 ROOT 的上级仓库，不在本 repo 内） |
| numpy | 2.5.3 |
| pandas | 3.0.6 |
| matplotlib | 3.11.2 |
| scipy | 1.18.1 |
| pymc / pytensor | 6.3.2 / 3.3.2 |
| arviz | 1.3.0 |
| nutpie | 0.16.11 |
| numba / llvmlite | 0.67.0 / 0.49.0 |

**smoke 验证实际使用**：`<ROOT上级>/.venv-nee-bayes/bin/python`（pandas 3.0.6 / matplotlib 3.11.2 /
numpy 2.5.3 / scipy 1.18.1），本 repo `bin/nee_recon.py` 在其中 8 张图全部重绘成功。

`requirements.lock.txt` 未包含 `python-docx`（DOCX 渲染脚本 `render_manuscript.py` 需要），也未包含
`pysam / skani / checkm2 / gtdbtk / iqtree` 等外部二进制。DOCX 渲染在本 repo 外层环境另装
`python-docx==1.2.0`（系统 python3 实测）。

## 外部二进制工具（Methods 声明版本；本 repo 不含，需在计算服务器环境复原）

| 工具 | 版本（Methods 声明） | 用途 | 本 repo 状态 |
|---|---|---|---|
| CheckM2 | v1.1.0 | 基因组质量（完整度/污染度） | 未随 repo；外部输入必需 |
| skani | v0.3.2 | ANI95 聚类 | 未随 repo；外部输入必需 |
| GTDB-Tk | release 226 | 分类与标记比对 | 未随 repo；外部输入必需 |
| IQ-TREE | v2.0.7 | 系统发育（LG+F+R10 / LG+F+I+G4） | 未随 repo；外部输入必需 |
| Prodigal | v2.6.3 | 蛋白预测 | 未随 repo |
| eggNOG-mapper | v2.1.15 | 功能注释 | 未随 repo |
| KOfamScan / HMMER | v1.3.0 / v3.3.2 | KOfam 注释 | 未随 repo |
| DIAMOND | (未标注版本) | 独立招募比对 | 未随 repo |
| geNomad | v1.9 | 移动元件上下文 | 未随 repo |
| MAFFT / HMMER / ESMFold / ProtNote / mDeepFRI v1.1.8 / DPFunc | 见 Methods | DPFQ008 候选表征 | 未随 repo；外部输入必需 |

## R 环境（部分主图/校正面板）

legacy 脚本中存在 R 版本（`207.make_ISME_figure2_lineage_mosaic.R`、`435–440` 校正重绘、
`255/256` DPFQ008 树分析）。R 版本**未在项目内冻结记录**（未联网核对）；`KO_panel_repair_20260926/R_sessionInfo.txt`
存在但未打包进本 repo。属 `external-input-required`。

## 数据来源层级

1. **公共测序数据**：权威 accession 见 `data/accessions/public_accessions.tsv`。
2. **大体积上游衍生输入**（注释矩阵、presence 表、矢量图）位于工作区外层的 `code/result_raw/`（**ROOT 之外**），
   以及 `NEE_revision_20260926/4x_*/host/` 各轮目录。**部分已打包**进 `data/upstream/`（见
   `data/upstream/SOURCES.tsv`，用相对来源标签 + sha256 记录，不含主机绝对路径）；其余大体积条目
   （FASTA/BAM/模型检查点/后验阵列）仍未打包，详见 `reproducibility_report.md`。
   本 repo 以相对来源标签 + sha256 记录 provenance，**不写任何主机绝对路径或主机名**。
