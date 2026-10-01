# 复现性报告（reproducibility_report.md）

**结论先行（诚实口径）**：本轮修的是**图件保真度**问题——上一版 CLI 曾把原始绘图逻辑
"重新发明"了一遍，导致 Fig2 只剩热图（丢失原始树与细菌/古菌分面）、Fig4 坐标轴错位且按位置
取列取到了字符串列、FigS3 把两种零模型混在一起且重复设标题、Fig3 丢失原始面板区分、
Fig5 panel d 丢掉了基因箭头。**现已全部改为直接移植作者端原始绘图代码**
（`scripts/figures/manuscript_figures.py`，文件头给出逐图来源索引），数值列一律**按列名**取用。

**十张图状态**：九张为 `regenerated-from-original-code`——**待 host 目视核对**；
Fig6 为 `verified`（原始矢量 + 52 轮已批准标签修正，SVG 与 52 轮归档**逐字节一致**）。
本仓库**不自我宣称**已通过视觉 QA；目视核对由 host 独立完成。

**仍未声称**：raw-to-paper 全流程复现。重型上游（拼接/注释/系统发育/贝叶斯/AI 模型）
依赖计算服务器与大体积输入，未随仓库打包，逐条点名见 D 节。

---

## A. 本轮实际执行并验证的项目

| 动作 | 命令（仓库根，全部相对路径） | 结果 |
|---|---|---|
| 图件输入完整性 | `python bin/nee_recon.py verify` | 0 缺失，退出码 0（Fig1..FigS4 共 10 行全 OK） |
| 十张图重绘 | `python bin/nee_recon.py plot all` | 10 png + 10 pdf + 1 svg，退出码 0，约 6 s |
| Fig6 标签修正链 | `python bin/nee_recon.py plot fig6` | `figure6_DPFQ008.svg` 与 52 轮归档 `figures/figure6_DPFQ008.svg` **sha256 相同** |
| 联合构型诊断 | `python bin/nee_recon.py joint` | 5 个输出与 `data/upstream/joint/reference_outputs/` **逐字节一致** |

环境：Python 3.12.12（`pandas 3.0.6 / matplotlib 3.11.2 / numpy 2.5.3`）；
Fig6 的 PDF/PNG 需要 `rsvg-convert`（脚本自动探测 PATH 与常见前缀，缺失回退 cairosvg）。

## B. 逐图状态（`manifest/figure_to_analysis.tsv` 为准）与移植来源

| 图 | 原始绘图代码来源 | 状态 |
|---|---|---|
| Fig1 | `r51/build_figures.py`（Figure 1 段） | regenerated-from-original-code |
| Fig2 | `r51/restore_tree_figure.py`（整文件：原始树 + 细菌/古菌分面 + 谱系标签） | regenerated-from-original-code |
| Fig3 | `r48/plot_diagnostics.py`（整文件，含原始统计注脚） | regenerated-from-original-code |
| Fig4 | `r51/build_figures.py`（Figure 4 段，列名取用） | regenerated-from-original-code |
| Fig5 | `r51/build_figures.py`（Figure 5 段，panel d 基因箭头） | regenerated-from-original-code |
| Fig6 | `r52/restore_candidate_figure.py` | verified（SVG byte-identical） |
| FigS1 | `r34/build_figures.py`（figure1_framework） | regenerated-from-original-code |
| FigS2 | `r34/build_figures.py`（figure3_bayesian_profiles） | regenerated-from-original-code |
| FigS3 | `r51/build_figures.py`（KO-null 段；仅 within-source 零模型 + 聚类计数） | regenerated-from-original-code |
| FigS4 | legacy `435.redraw_corrected_WGS_figure_20260927.R`（(A\|B)/C 版式等价比移植） | regenerated-from-original-code |

**状态一致性**：`regenerated-from-original-code` 9 / `verified` 1 / `external-input-required` 0；
视觉 QA 状态：**待 host 目视核对**（本仓库不作视觉 QA 自证）。

## C. 本轮修复的具体保真度缺陷

1. **Fig2**：旧实现只有一张热图；现按 `restore_tree_figure.py` 恢复**双侧原始树**
   （`data/derived/{bacteria,archaea}_tree_nodes|edges.tsv`）+ 谱系标签列 + 细菌/古菌分面，
   色标 0–1 与稿件一致。
2. **Fig4**：旧实现 `cod.columns` 位置切分导致来源/系统轴对调，且 `ep.iloc[:,1/2]` 取到了
   `cohort/rank` 等**字符串元数据**当作数值轴（出现巨大 x 轴标签）。现改为
   `aa.source_or_site_partial_r2` / `aa.df_adjusted_partial_r2` **按列名**取用，面板方向与稿一致。
3. **FigS3**：旧实现把 `global` 与 `within-source` 两种零模型**混在一个散点里**且重复设标题；
   现按原代码只画 `within_source_prevalence_preserving`，a/b 分别为细菌/古菌，
   c 面板同时给出 global 与 within-source 的 BH 显著计数（含分母标签与图例）。
4. **Fig3**：旧实现只用衍生表近似；现整文件移植 `r48/plot_diagnostics.py`，
   保留 A/B/C 三面板区分（4 队列 × 2 权重方案、7 个模块、3 行面板敏感性）与原始统计注脚。
5. **Fig5**：旧实现 panel d 用方块点代替基因；现恢复原始 `FancyArrow` 箭头、基因名
   （kdpA/kdpB/kdpC、proV/proW/proX、ectA/ectB/ectC、betA/betB）、配色与 genome·source 注释。
6. **FigS4**：旧实现是单行三面板，与原始 R 版式不符；现按 `(A|B)/C` 版式重绘
   （A 策略检出、B 跨库计数、C 门槛敏感），标签/配色/形状映射（15/16/17/18）与原 R 一致，
   并保留三条原始断言（`represented_KOs` 和 = 65、跨库计数一致、12 行门槛表）。

## D. 仍未打包的输入（点名，未伪造）

1. **服务器注释/系统发育管线**（CheckM2 / Prodigal-meta / Kofam / GTDB-Tk r226 / IQ-TREE）：
   已提供去私有化参数化公开版 `scripts/server/`，标注 integration-not-retested。
2. **PyMC 后验阵列**（`*.nc`）：FigS2 使用原始 34 脚本内固定的 5 行归档值。
3. **DPFQ008 AI 预测输入**：ProtNote seed-42 检查点、mDeepFRI v1.1.8、DPFunc 适配输入
   （ESM2-650M / ESMFold / 22,369 维 InterPro 向量）——仅保留分数表 S3 与常规证据表。
4. **原始 reads / DIAMOND recruitment 比对**（FigS4 上游）：仅保留逐 KO/逐库汇总表。
5. `CITATION.cff` 的 `doi` 与 `license`：作者确认前**刻意留空**，未编造。

## E. 边界与保留的原始口径

- 本 repo **不产出** DOCX/PDF 稿件；`render_manuscript.py` 保留但需 `python-docx` 与稿件 `.md` 输入。
- **保留原始 caveat 文本**：Fig3 的 "Conditional on observed catalogues and coarse lineage blocks /
  equal subgroup weighting is a post-result diagnostic" 统计注脚、联合诊断合同的
  "Cox/Cut co-detection is NOT a confirmed CO-oxidation marker"、以及
  `data/upstream/fig6_dpfq008/INPUT_SOURCES.md` 的 "method-specific scores that are not comparable
  to one another"（AI 概率不可比）——均在原处保留，未删改。
- **未声称** "full raw-to-paper reproduction"；本 repo 为开发性/复现性代码发布，
  不代表已有正式发表文章。
