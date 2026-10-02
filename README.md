# Mars-analogue genomic strategies

Code and small source tables for the manuscript *Lineage structure and molecular diversity of microbial strategies across terrestrial Mars analogues*. This initial research-code release covers six main figures, five supplementary figures, and the frozen joint-configuration diagnostic. It does not reproduce the complete raw-read-to-paper workflow: large sequence inputs, reference databases and model weights are external. Detailed limitations and provenance are recorded below and in `reproducibility_report.md`.

```bash
python bin/nee_recon.py verify
python bin/nee_recon.py plot all
python bin/nee_recon.py joint
```

The added-cohort workflow is `scripts/server/run_added_cohorts.py`; configure external tools and databases using the adjacent example JSON. It processes a frozen 607-record intake only after the current download has completed and all official checksums pass. These records are **candidates**, not accepted genomes or independent ecological samples. ANI and source/site deduplication are subsequent steps. The three corrected RefSeq URLs and complete download manifest are included under `data/accessions/`.

## Expanded-cohort phylogeny

The [round 60 workflow](experiments/final_phylogeny/README.md) provides the frozen 4,477-genome input manifest and quality-specific tree inference code. Marker extraction has started; completed expanded trees and revised ecological conclusions are not yet available. External sequences and reference databases are required.

## Resource-niche feasibility pilot

The optional [resource-niche annotation pilot](experiments/resource_niche_pilot/README.md) reproduces a 16-MAG parser-sensitivity experiment and the Mackay missing-measurement audit. These are exploratory method checks, outside the manuscript results; they do not test competition or environmental mechanisms. Pinned upstream files are fetched and hash-checked separately.

## Earlier ten-figure reproduction notes

The detailed notes below describe the inherited v6 figures. The v7 addition, FigS5, is generated directly from the frozen four-pair joint-state analysis. Current availability is determined by the figure manifest, not the historical ten-figure count.

# 论文代码与复现包 —— v6（manuscript_v6, 6 主图 + 4 补图）

Mars analogue 微生物策略论文（Ruxin Sun, Shaocheng Yan ... Yuming Fu）的**公开可交付代码仓库**。
本仓库把权威稿 `manuscript_v6.md` 的 **6 张主图 + 4 张补充图**及其核心结论，逐条映射到**真实存在**
的脚本与**已审计的小体积**数据（合计 `data/upstream` 约 12 MB，逐项见 `data/upstream/SOURCES.tsv`）。

> **重要边界（请先读）**：本 repo 是**图件与核心表的复现层**，**不是** "raw-to-paper" 全流程复现。
> 原始 MAG 拼接、注释（Prodigal/eggNOG/KOfam）、系统发育（GTDB-Tk/IQ-TREE）、贝叶斯拟合（PyMC）
> 等重型上游步骤依赖计算服务器与大体积输入（FASTA/BAM/模型检查点/后验阵列），**未随本仓库打包**。
> 本仓库内**十张图全部可由 CLI 重新产出**，且绘图**不是重新发明**：每一步都取自作者端
> **原始绘图代码**（`scripts/figures/manuscript_figures.py` 头部给出逐图来源），数值列一律按列名取用。
> 其中 Fig6 采用已批准的「原始矢量 + 标签修正」链（SVG 与 52 轮归档逐字节一致）、
> FigS4 按原 R 脚本的 (A|B)/C 版式从逐 KO/逐库源表重绘。每张图/每条结论的状态逐条标注，见
> `manifest/figure_to_analysis.tsv` 与 `manifest/claims_to_evidence.tsv`。**目视核对由 host 独立完成**，
> 本仓库不自我宣称已通过视觉 QA。
>
> 本仓库为**开发性/复现性代码发布**，不代表已有正式发表文章。

## 1. 目录结构

```
code_repo/
├── README.md                     # 本文件
├── CITATION.cff                  # 引用元数据（无 DOI/license，见 docs/）
├── .gitignore                    # 忽略可再生成产物与敏感/大体积内容
├── bin/
│   └── nee_recon.py              # 可复用 CLI（相对路径，--input-root/--output-root）
├── manifest/
│   ├── figure_to_analysis.tsv    # ★ 十张图 → 数据/脚本/命令/输出/状态
│   ├── claims_to_evidence.tsv    # ★ 25 条核心结论 → 值/证据文件/状态
│   └── provenance_hashes.tsv     # 本 repo 每个文件的 sha256 + 相对来源标签
├── data/
│   ├── derived/                  # 小体积图件数值输入（23 个 TSV）
│   ├── accessions/               # 公共 accession / 站点 / 候选证据元数据
│   └── upstream/                 # ★ 已审计小体积上游输入（约 12 MB）
│       ├── SOURCES.tsv           #   每个文件的相对来源标签 + sha256（无主机路径）
│       ├── panel_repair/         #   KO 面板修正如 amended portfolio/recurrence/矢量树源
│       ├── cohorts/              #   49/48 轮队列状态与代表元数据（已剥离私有路径列）
│       ├── context/              #   51 轮靶蛋白/坐标/紧凑位点表
│       ├── fig6_dpfq008/         #   Fig6 原始矢量、源树、匹配零模型 + INPUT_SOURCES.md
│       ├── figs4/                #   FigS4 逐 KO/逐库/阈值源表 + WGS 证据
│       └── joint/                #   联合构型诊断契约、输入说明与参考输出
├── scripts/
│   ├── figures/                  # ★ 绘图实现层（十张图的唯一绘图代码）
│   │   └── manuscript_figures.py #   原始绘图代码移植版；文件头给出逐图来源索引
│   ├── joint/                    # ★ 可移植联合构型支持 CLI
│   ├── server/                   # ★ 去私有化的服务器脚本（参数化 --workdir/--dbroot）
│   └── pipeline/                 # 48–52、15、34 轮真实管线脚本（原样复制）
│       ├── r15/ r34/ r48/ r49/ r50/ r51/ r52/
├── legacy/                       # 上游原始脚本（provenance 保留，不可直接运行）
│   ├── dpfq008/ primary_figures/ corrected_panels/
├── env/
│   ├── requirements.lock.txt     # 真实 pip freeze（15 轮）
│   └── environment_notes.md      # Python/R/二进制工具版本与缺口
├── docs/citation_metadata.md     # 标题/作者/资助（照录稿件）
└── reproducibility_report.md     # ★ 逐条缺失项点名与修复状态
```

## 2. 十张图一览

| 图 | 标题 | 状态 |
|---|---|---|
| Fig1 | Catalogue composition and encoded strategy profiles | regenerated-from-original-code（待 host 目视核对） |
| Fig2 | Encoded strategies form heterogeneous mosaics across lineages | regenerated-from-original-code（含原始树 + 细菌/古菌分面；待 host 目视核对） |
| Fig3 | The maintenance–energy source contrast depends on uncertainty and aggregation | regenerated-from-original-code（A/B/C 三面板 + 原始统计注脚；待 host 目视核对） |
| Fig4 | Specific molecular combinations carry unequal source associations | regenerated-from-original-code（数值列按名取用；待 host 目视核对） |
| Fig5 | Local genomic organization supports some combinations | regenerated-from-original-code（panel d 原始基因箭头；待 host 目视核对） |
| Fig6 | AI-assisted prioritization of DPFQ008 | verified（原始矢量标签修正；SVG 与 52 轮输出逐字节一致） |
| FigS1 | Constraint–resource–strategy framework (conceptual) | regenerated-from-original-code（待 host 目视核对） |
| FigS2 | Existing Bayesian source-profile contrast | regenerated-from-original-code（待 host 目视核对） |
| FigS3 | Per-KO source-preserving null results | regenerated-from-original-code（仅 within-source 零模型 + 聚类计数面板；待 host 目视核对） |
| FigS4 | Independent sequence occurrence | regenerated-from-original-code（按原 R 脚本 (A\|B)/C 版式重绘；待 host 目视核对） |

逐图的**数据源、脚本、命令、输出、状态**见 `manifest/figure_to_analysis.tsv`。
**状态口径**：九张为「由原始绘图代码重绘、待 host 目视核对」，Fig6 为「原始矢量标签修正」。
本仓库**不自我宣称**已通过视觉 QA——目视核对由 host 独立完成。

## 3. 如何运行

### 3.1 环境
需要一个含 `pandas / matplotlib / numpy` 的 Python（本项目实测：pandas 3.0.6,
matplotlib 3.11.2, numpy 2.5.3 — 见 `env/requirements.lock.txt`）。Fig6 的 PDF/PNG 输出额外需要
`rsvg-convert`（脚本自动探测 `PATH` 与常见安装前缀；缺失时回退 `cairosvg`，再缺失则仅保留 SVG）。

```bash
python3 -m venv .venv && . .venv/bin/activate && pip install -r env/requirements.lock.txt
```

### 3.2 CLI（全部相对路径）

```bash
# 列出十张图与状态
python bin/nee_recon.py list

# 校验已打包图件输入是否存在
python bin/nee_recon.py verify

# 重绘单张 / 全部十张图（默认输出到 <repo>/figures_out）
python bin/nee_recon.py plot fig6
python bin/nee_recon.py plot all

# 联合构型支持诊断（可移植；写入 <out>/joint）
python bin/nee_recon.py joint

# 自定义路径（发布/CI 推荐）
python bin/nee_recon.py --input-root /path/to/repo --output-root /tmp/out plot fig4
```

`plot` 调用 `scripts/figures/manuscript_figures.py`（**原始绘图代码移植版**）重绘全部十张图
（`fig1..fig6, figS1..figS4`），输入仅为本仓库已打包的小体积衍生表 / 源表 / 原始矢量。
`joint` 运行 `scripts/joint/joint_configuration_support.py`，其输出与作者端
`data/upstream/joint/reference_outputs/` 的参考输出**逐字节一致**。

### 3.3 真实管线脚本（`scripts/pipeline/`、`scripts/server/`）

`scripts/pipeline/` 是 48–52 / 15 / 34 轮的**原始脚本**（原样复制），含硬编码的
`Path(__file__)` 相对定位，指向各轮 `host/`，**需要对应的上游输入**才能整体运行。
`scripts/server/` 提供**去私有化、参数化**的公开版本（`--workdir/--dbroot/--gtdb-*` 等），
供在自有服务器上复原；这些脚本**未在本 repo 端到端重跑**（标注 integration-not-retested），
只作为归档分析代码。

## 4. 数据与 provenance

- `data/derived/`：各图件的数值输入表（均 < 40 KB，合计约 250 KB）。
- `data/upstream/`：已审计的小体积上游输入（约 12 MB）；`data/upstream/SOURCES.tsv` 逐文件给出
  **相对来源标签 + sha256**。
- `data/accessions/`：公共 accession、站点范围、候选证据（公共元数据，小体积）。
- `manifest/provenance_hashes.tsv`：本 repo 每个文件的 sha256 + 相对来源标签；`legacy/` 保留原脚本用于溯源。
- **隐私**：本 repo 不写任何主机绝对路径或主机名；`data/upstream/cohorts/` 中的新队列状态/代表表
  已**剥离两张私有路径列**（MAG FASTA 路径与来源清单），保留全部生物 id 列。

## 5. 复现状态与缺口

- 十图状态：九张 `regenerated-from-original-code`（**待 host 目视核对**）+ Fig6 `verified`
  （原始矢量标签修正，SVG 与 52 轮归档逐字节一致）。状态口径见 `manifest/figure_to_analysis.tsv`。
- 绘图实现集中在 `scripts/figures/manuscript_figures.py`，与稿件版面/配色/面板划分/统计注脚一致；
  数值一律按列名取用，不使用位置索引。
- `data/upstream/joint/` 的联合构型诊断**可移植重跑并通过**（5 个输出与参考逐字节一致）；
  判据（4 对 Cox/Cut 共检=**非**确认 CO、门槛 2/3/5、group≥3）与原始合同未改动。
- 仍未打包（点名见 `reproducibility_report.md`）：服务器注释/系统发育管线、PyMC 后验阵列、
  DPFQ008 AI 模型检查点、原始 reads/recruitment 比对。
- **未声称** "full raw-to-paper reproduction"；本 repo 不包含稿件、PDF、FASTA/BAM 或任何含 token 的日志/私有配置。
