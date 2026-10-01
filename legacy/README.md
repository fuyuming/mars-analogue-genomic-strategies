# legacy/ —— 上游原始脚本（provenance 保留）

本目录保存被 48–52 轮交付管线**引用**、但不在最近轮次目录内的原始脚本，仅用于溯源：
它们产生 51/52 轮图件所依赖的中间产物或原始矢量。**未做任何修改**，逐文件 sha256 见
`../manifest/provenance_hashes.tsv`（`repo_path` 以 `legacy/` 开头）。

来源根：`<workspace>/code`（**在 ROOT 之外**，编号脚本 1–442）。这些脚本多为 R / Python / shell，
依赖计算服务器与 `code/result_raw/` 下的大体积输入，**在本 repo 内不可直接运行**。
（本 repo 用相对来源标签 + sha256 记录 provenance，不写任何主机绝对路径。）

## legacy/dpfq008/ —— DPFQ008 候选家族分析（图 6 上游）
| 文件 | 作用 | 对应结论 |
|---|---|---|
| 250.extract_DPFQ008_family_architecture.py | 提取代表序列结构（CXXCH、InterPro 区间） | C22 |
| 251.test_DPFQ008_matched_context_enrichment.py | within-MAG 匹配零模型（100,000 次置换） | C24 |
| 253.extract_DPFQ008_hmm_hits.py | 严格 profile-HMM 命中（63 蛋白/61 MAG/58 ANI95） | C22 |
| 255.analyze_DPFQ008_gene_host_tree_congruence.R | 基因树 vs 宿主树距离相关（Mantel） | C23 |
| 256.make_ISME_figure3_DPFQ008_lineage_context.R | 生成图 6 原始矢量 SVG | 图 6 |
| 257.run_DPFQ008_hmm_seed_cv_hpc.sh | 5 折种子敏感性（原文件名含主机名，已中性化） | C22 |
| 260.make_DPFQ008_validation_supplement_table.py | 生成 Supplementary Table S3 | 表 S3 |

## legacy/primary_figures/ —— 早期主图脚本
`207`（lineage mosaic，图 2 上游）、`209`（strict tree 集成）、`216`（mobile context）、
`224`（AI neighborhood）、`25`/`38`（Nature/ISME 4213 版主图）。

## legacy/corrected_panels/ —— KO 面板修正与校正重绘（图 S4、图 2 校正）
`411/413`（core65 / restored74 面板重建）、`435`（WGS 图 = 图 S4）、`436`（lineage 图）、
`437`（portfolio supplement）、`439`（source 图）、`440`（mobile 面板导出）。

> 说明：这些脚本引用的输入（如 `code/result_raw/KO_panel_repair_20260926/...`、
> `code/result_raw/ISME_Figure3_DPFQ008_lineage_context_20260729_v11/...`）**未随本 repo 打包**。
> 缺哪些、差什么见 `../reproducibility_report.md`。
