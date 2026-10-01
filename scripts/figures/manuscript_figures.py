# -*- coding: utf-8 -*-
"""manuscript_figures —— 十张图的**忠实重绘**（原始绘图代码移植版）。

本模块不是新画的图，而是把作者端**原始绘图代码段**逐段移植为可移植、相对路径的版本，
输出与稿件 `manuscript_v6/assets/` 中的图件一一对应（6 主图 + 4 补图）。

来源索引（function -> original code）
  fig1_cohort_strategies       <- 51_evidence_figures_20261001/build_figures.py  (Figure 1 段)
  fig2_lineage_mosaic          <- 51_evidence_figures_20261001/restore_tree_figure.py (整文件)
  fig3_portfolio_uncertainty   <- 48_MAG_portfolio_inference_20261001/plot_diagnostics.py (整文件)
  fig4_molecular_configurations<- 51_evidence_figures_20261001/build_figures.py  (Figure 4 段)
  fig5_dryland_context         <- 51_evidence_figures_20261001/build_figures.py  (Figure 5 段)
  figs1_framework              <- 34_integrated_manuscript_20260929/build_figures.py (figure1_framework)
  figs2_bayesian_profiles      <- 34_integrated_manuscript_20260929/build_figures.py (figure3_bayesian_profiles)
  figs3_null_details           <- 51_evidence_figures_20261001/build_figures.py  (Figure S3 / KO-null 段)
  figs4_independent_occurrence <- legacy/corrected_panels/435.redraw_corrected_WGS_figure_20260927.R
  fig6_dpfq008                 <- 52_framework_integration_20261001/restore_candidate_figure.py

规则
- 数值列一律**按列名**访问（不用位置索引）。
- 全部使用相对路径：`root` = 仓库根，输入取自 data/derived 与 data/upstream。
- 不改动任何科学阈值、不重跑重型分析；绘图参数（配色/面板布局/标签/统计注释）与原代码一致。
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, Patch, FancyBboxPatch

TEAL = "#176B80"
ORANGE = "#BB6C32"
GREY = "#B3B9BC"

SOURCES = ["qaidam_basin_mmag_1773", "alaska_permafrost_reference", "stordalen_mire_2019_hybrid_mags",
           "atacama_salt_crust_prjna351262", "atacama_halite_rainfall_prjna484015",
           "mauna_loa_lava_tube_fishman_2023", "australian_basalt_lava_tubes_bay_2025"]
LABELS = ["Qaidam", "Alaska", "Stordalen", "Atacama crust", "Atacama halite", "Mauna Loa", "Australian caves"]
SYSTEMS = ["KdpABC", "KdpDE", "ProVWX", "EctABC", "EctD", "BetAB"]
MODULE_NAMES = {"cold_protein_quality": "Cold / protein quality",
                "dna_repair_radiation": "DNA repair / protection",
                "dormancy_resuscitation": "Dormancy / stringent",
                "osmotic_desiccation_salt": "Solute / ion homeostasis",
                "oxidative_redox": "Oxidative / redox",
                "light_energy": "Light-associated",
                "sulfur_chemolithotrophy": "Sulfur metabolism",
                "trace_gas_energy": "H2 / CO / carbon markers",
                "biofilm_eps_surface": "Surface / EPS"}


def _rc(**over):
    p = {"font.family": "Arial", "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
         "xtick.labelsize": 10, "ytick.labelsize": 10, "pdf.fonttype": 42, "svg.fonttype": "none",
         "axes.spines.top": False, "axes.spines.right": False}
    p.update(over)
    plt.rcParams.update(p)
    return plt


def _save(fig, out, name, exts=("png", "pdf"), dpi=230, bbox=True):
    out.mkdir(parents=True, exist_ok=True)
    made = []
    kw = {"dpi": dpi, "facecolor": "white"}
    if bbox:
        kw["bbox_inches"] = "tight"
    for ext in exts:
        p = out / f"{name}.{ext}"
        fig.savefig(p, **kw)
        made.append(str(p))
    plt.close(fig)
    return made


def _tsv(p):
    return pd.read_csv(p, sep="\t")


def _title(ax, letter, s, pad=12):
    ax.set_title(letter + "  " + s, loc="left", fontweight="bold", pad=pad)


# --------------------------------------------------------------------------------------
# Figure 1 —— 51/build_figures.py Figure 1 段（逐段移植）
# --------------------------------------------------------------------------------------
def fig1_cohort_strategies(root, out):
    plt = _rc()
    U = root / "data/upstream"
    b = pd.read_csv(U / "panel_repair/amended_inputs/core65_MAG_module_breadth.tsv.gz", sep="\t")
    meta = b.drop_duplicates("genome_id")
    rep = meta[meta["species_representative"]].copy()
    counts = meta.groupby("dataset_id").agg(primary=("genome_id", "size"),
                                            strict=("passes_strict_90_5", "sum")).reindex(SOURCES)
    counts["ANI95"] = rep.groupby("dataset_id").size()
    top = list(rep.phylum.value_counts().head(5).index) + ["Halobacteriota"]
    rep["display_phylum"] = np.where(rep.phylum.isin(top), rep.phylum, "Other / unclassified")
    ct = pd.crosstab(rep.dataset_id, rep.display_phylum).reindex(SOURCES)
    ct = ct.reindex(columns=list(top) + ["Other / unclassified"])
    frac = ct.div(ct.sum(axis=1), axis=0)
    co = _tsv(U / "cohorts/core65_primary_cohort.tsv")
    means = (b[b.genome_id.isin(co.genome_id)].groupby(["dataset_id", "module"])
             .module_breadth.mean().unstack().reindex(index=SOURCES, columns=list(MODULE_NAMES)))

    fig = plt.figure(figsize=(10, 7.5), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
    a = fig.add_subplot(gs[0, 0]); c = fig.add_subplot(gs[0, 1]); e = fig.add_subplot(gs[1, :])
    y = np.arange(7)
    a.barh(y, counts.primary, color="#DCE3E5", label="Primary MAGs")
    a.barh(y, counts.strict, color=TEAL, label="Strict subset")
    a.plot(counts.ANI95, y, "D", ms=4.5, color=ORANGE, label="ANI95 representatives")
    for i, v in enumerate(counts.primary):
        a.text(v + 20, i, str(v), va="center", fontsize=9)
    a.set(yticks=y, yticklabels=LABELS, xlabel="Genomes", xlim=(0, 1900))
    a.invert_yaxis()
    a.legend(frameon=False, fontsize=9, loc="center right", bbox_to_anchor=(1, .40))
    _title(a, "a", "Catalogue representation")
    left = np.zeros(7)
    colors = ["#176B80", "#5F9F9D", "#DEA158", "#A67EAB", "#7C8C53", "#887455", "#D9D9D9"]
    for col, colr in zip(frac.columns, colors):
        c.barh(y, frac[col], left=left, color=colr, label=col)
        left += frac[col].values
    c.set(yticks=y, yticklabels=[], xlabel="Fraction of ANI95 representatives", xlim=(0, 1))
    c.invert_yaxis()
    _title(c, "b", "Taxonomic composition")
    c.legend(frameon=False, fontsize=9, bbox_to_anchor=(1.01, .98), loc="upper left")
    im = e.imshow(means.values.T, aspect="auto", vmin=0, vmax=1, cmap="YlGnBu")
    e.set(xticks=y, xticklabels=LABELS, yticks=range(9),
          yticklabels=[MODULE_NAMES[m] for m in means.columns])
    e.tick_params(axis="x", labelrotation=18)
    for i in range(9):
        for j in range(7):
            e.text(j, i, ("<0.01" if 0 < means.iloc[j, i] < .005 else f"{means.iloc[j, i]:.2f}"),
                   ha="center", va="center", fontsize=9, color="white" if means.iloc[j, i] > .6 else "black")
    _title(e, "c", "Encoded module breadth in the bacterial comparison")
    fig.colorbar(im, ax=e, shrink=.8, label="Mean fraction of specified KOs")
    return _save(fig, out, "figure1_cohort_strategies")


# --------------------------------------------------------------------------------------
# Figure 2 —— 51/restore_tree_figure.py（整文件移植：原始树 + 细菌/古菌分面 + 谱系标签）
# --------------------------------------------------------------------------------------
def fig2_lineage_mosaic(root, out):
    _rc()
    U = root / "data/upstream"
    D = root / "data/derived"
    S = U / "panel_repair/figure2_source"
    means = pd.read_csv(S / "Figure2B_lineage_module_breadth.csv")
    order = pd.read_csv(S / "Figure2B_tree_tip_order.csv")
    meta = pd.read_csv(S / "monophyletic_lineage_summary.csv").set_index("lineage_id")
    modules = ["cold_protein_quality", "dna_repair_radiation", "dormancy_resuscitation",
               "osmotic_desiccation_salt", "oxidative_redox", "light_energy",
               "sulfur_chemolithotrophy", "trace_gas_energy", "biofilm_eps_surface"]
    labels = ["Cold / protein\nquality", "DNA repair /\nprotection", "Dormancy /\nstringent",
              "Solute / ion\nhomeostasis", "Oxidative /\nredox", "Light-\nassociated",
              "Sulfur\nmetabolism", "H2 / CO /\ncarbon", "Surface /\nEPS"]
    fig = plt.figure(figsize=(10, 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 3, height_ratios=[4.6, 1.2], width_ratios=[1.3, 1.6, 5.2],
                          hspace=.14, wspace=.025)
    for i, dom in enumerate(["Bacteria", "Archaea"]):
        ax = fig.add_subplot(gs[i, 0]); la = fig.add_subplot(gs[i, 1]); hm = fig.add_subplot(gs[i, 2])
        ids = order[order.domain == dom].sort_values("y", ascending=False).lineage_id.tolist()
        pos = {k: j for j, k in enumerate(ids)}
        nodes = pd.read_csv(D / f"{dom.lower()}_tree_nodes.tsv", sep="\t").fillna("")
        edges = pd.read_csv(D / f"{dom.lower()}_tree_edges.tsv", sep="\t")
        depth = nodes.set_index("node").depth.to_dict()
        nn = nodes.set_index("node").label.to_dict()
        child = edges.groupby("parent").child.apply(list).to_dict()
        ys = {}

        def calc(n):
            if n in ys:
                return ys[n]
            ys[n] = pos[nn[n]] if nn[n] else float(np.mean([calc(c) for c in child[n]]))
            return ys[n]

        for n in nodes.node:
            calc(n)
        for p, cs in child.items():
            ax.plot([depth[p]] * 2, [min(ys[c] for c in cs), max(ys[c] for c in cs)], color="#69777E", lw=.6)
            for c in cs:
                ax.plot([depth[p], depth[c]], [ys[c]] * 2, color="#69777E", lw=.6)
        ax.set(ylim=(len(ids) - .5, -.5), xlabel="Branch length", yticks=[])
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="x", labelsize=9)
        ax.set_title(("a" if i == 0 else "b") + "  " + dom + "\n" + str(len(ids)) + " lineage blocks",
                     loc="left", fontsize=13, fontweight="bold")
        la.set(xlim=(0, 1), ylim=(len(ids) - .5, -.5)); la.axis("off")
        chosen = []
        if i == 0:
            for k in meta.loc[ids].sort_values("genomes", ascending=False).index:
                if meta.loc[k, "genomes"] >= 20 and all(abs(pos[k] - pos[x]) >= 4 for x in chosen):
                    chosen.append(k)
        else:
            chosen = ids
        for k in chosen:
            la.text(.98, pos[k], str(meta.loc[k, "collapsed_name"]) + " (" + str(meta.loc[k, "genomes"]) + ")",
                    ha="right", va="center", fontsize=10.5)
        m = means[means.domain == dom].pivot(index="lineage_id", columns="module",
                                             values="mean_module_breadth").loc[ids, modules]
        im = hm.imshow(m, aspect="auto", vmin=0, vmax=1, cmap="YlGnBu", interpolation="nearest")
        hm.set(yticks=[], xticks=range(9))
        hm.set_xticklabels(labels if i == 1 else [], rotation=45, ha="right")
        hm.tick_params(axis="both", length=0)
        hm.spines[:].set_visible(False)
        if i == 0:
            hm.set_title("Mean encoded marker breadth by lineage", loc="left", fontsize=13, fontweight="bold")
        for x in [4.5, 7.5]:
            hm.axvline(x, color="white", lw=2)
    fig.colorbar(im, ax=fig.axes, location="right", shrink=.65, pad=.02, label="Mean fraction of specified KOs")
    return _save(fig, out, "figure2_lineage_mosaic")


# --------------------------------------------------------------------------------------
# Figure 3 —— 48/plot_diagnostics.py（整文件移植；保留原始统计注脚）
# --------------------------------------------------------------------------------------
def fig3_portfolio_uncertainty(root, out):
    plt = _rc(**{"font.size": 9})
    D = root / "data/derived"
    d = _tsv(D / "paired_uncertainty.tsv")
    m = _tsv(D / "module_decomposition.tsv")
    eq = _tsv(D / "equal_submodule_uncertainty.tsv")
    fig, ax = plt.subplots(1, 3, figsize=(13.7, 4.6), gridspec_kw={"width_ratios": [1.05, 1.1, 1.15]})
    labels = []
    for i, (branch, strict) in enumerate([("core65", False), ("core65", True),
                                          ("restored74", False), ("restored74", True)]):
        labels.append(branch + (" strict" if strict else " primary"))
        for scheme, off, color in [("lineage_BB", .12, "#1b6a82"),
                                   ("lineage_source_product_sensitivity", -.12, "#b26b36")]:
            z = d[(d.branch == branch) & (d.strict == strict) & (d.scheme == scheme)
                  & (d.comparison == "full")].iloc[0]
            ax[0].plot([z.delta_q025 * 100, z.delta_q975 * 100], [i + off] * 2, color=color, lw=2)
            ax[0].plot(z.delta_median * 100, i + off, "o", color=color, ms=4,
                       label=("Lineage weights" if off > 0 else "Lineage × source weights") if i == 0 else None)
    ax[0].set(yticks=range(4), yticklabels=labels, xlabel="Maintenance − energy partial R² (pp)",
              title="A  Paired uncertainty")
    ax[0].set_ylim(4.3, -.5); ax[0].legend(frameon=False, fontsize=8, loc="lower right")
    ax[0].axvline(0, c=".6", lw=.8, ls="--")
    z = m[(m.branch == "core65") & ~m.strict]
    names = ["Cold-associated", "DNA repair/protection", "Sporulation/stringent", "Solute/ion homeostasis",
             "Oxidative/redox", "Sulfur metabolism", "H₂/CO/carbon markers"]
    ax[1].barh(range(7), 100 * z.source_partial_r2, color=["#1b6a82"] * 5 + ["#b26b36"] * 2)
    ax[1].set(yticks=range(7), yticklabels=names, xlabel="Module source partial R² (%)",
              title="B  Module decomposition")
    ax[1].invert_yaxis()
    for i, (comp, label) in enumerate([("full", "Original module weighting"),
                                       ("drop_osmotic_desiccation_salt", "Omit solute/ion module"),
                                       ("equal_osmotic_submodules", "Equal solute/ion subgroup weights")]):
        z = (eq if i == 2 else d)
        z = z[(z.branch == "core65") & ~z.strict & (z.scheme == "lineage_BB") & (z.comparison == comp)].iloc[0]
        ax[2].plot([z.delta_q025 * 100, z.delta_q975 * 100], [i] * 2, c="#1b6a82", lw=2)
        ax[2].plot(z.delta_median * 100, i, "o", c="#1b6a82")
        ax[2].text(z.delta_q025 * 100, i - .14, label, fontsize=8)
    ax[2].set(yticks=[], xlabel="Maintenance − energy partial R² (pp)", title="C  Panel sensitivity (core65)")
    ax[2].set_ylim(2.6, -.6); ax[2].axvline(0, c=".6", ls="--", lw=.8)
    fig.text(.01, .025, "Intervals: 95% reweighting ranges; points: medians. Conditional on observed catalogues "
                        "and coarse lineage blocks.\nEqual subgroup weighting is a post-result diagnostic, "
                        "not a replacement primary endpoint.", fontsize=8)
    fig.tight_layout(rect=[0, .105, 1, 1], w_pad=2)
    return _save(fig, out, "figure3_portfolio_uncertainty", dpi=200, bbox=False)


# --------------------------------------------------------------------------------------
# Figure 4 —— 51/build_figures.py Figure 4 段（数值列按名访问）
# --------------------------------------------------------------------------------------
def fig4_molecular_configurations(root, out):
    plt = _rc()
    U = root / "data/upstream"
    D = root / "data/derived"
    old = _tsv(U / "cohorts/original_primary_system_states.tsv")
    m = old.groupby("dataset_id")[SYSTEMS].agg(lambda x: (x == "panel_complete").mean()).reindex(SOURCES)
    fig = plt.figure(figsize=(10, 6.3), layout="constrained")
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1.2, 1])
    a = fig.add_subplot(gs[0, :]); c = fig.add_subplot(gs[1, 0]); e = fig.add_subplot(gs[1, 1])
    im = a.imshow(m.values.T * 100, aspect="auto", vmin=0, vmax=100, cmap="YlGnBu")
    a.set(xticks=range(7),
          xticklabels=[label + "\n(n=" + str(int((old.dataset_id == src).sum())) + ")"
                       for src, label in zip(SOURCES, LABELS)],
          yticks=range(6), yticklabels=SYSTEMS)
    for i in range(6):
        for j in range(7):
            a.text(j, i, f"{m.iloc[j, i] * 100:.1f}", ha="center", va="center",
                   color="white" if m.iloc[j, i] > .6 else "black", fontsize=10)
    _title(a, "a", "Specific marker combinations vary across original sources")
    fig.colorbar(im, ax=a, shrink=.9, label="Genomes with all specified markers (%)")
    assoc = _tsv(U / "cohorts/lineage_group_associations.tsv")
    aa = (assoc[(assoc.cohort == "original_primary") & (assoc["rank"] == "family")]
          .set_index("endpoint").loc[[s + "_codetected" for s in SYSTEMS]])
    c.plot(aa.source_or_site_partial_r2 * 100, range(6), "o", color=TEAL, label="Raw partial R²")
    c.plot(aa.df_adjusted_partial_r2 * 100, range(6), "x", color=ORANGE, label="Degrees-of-freedom adjusted")
    c.axvline(0, color=".75", lw=.7)
    c.set(yticks=range(6), yticklabels=SYSTEMS, xlabel="Source partial R² (%)", xlim=(-2, 25))
    c.invert_yaxis(); _title(c, "b", "Within 30 supported families")
    c.legend(frameon=False, fontsize=9, loc="lower right")
    cx = _tsv(D / "system_context_summary.tsv")
    cc = (cx[(cx.cohort == "original_strict") & (cx.metric == "order_compact")]
          .set_index("system").loc[[s for s in SYSTEMS if s != "EctD"]])
    e.barh(range(5), cc.fraction, color=TEAL, height=.55)
    for i, row in enumerate(cc.itertuples()):
        e.text(row.fraction + .025, i, f"{row.n}/{row.denominator}", va="center", fontsize=10)
    e.set(yticks=range(5), yticklabels=cc.index, xlabel="Fraction among co-detected genomes", xlim=(0, 1.18))
    e.invert_yaxis(); _title(e, "c", "Original CDS-order proximity proxy")
    return _save(fig, out, "figure4_molecular_configurations")


# --------------------------------------------------------------------------------------
# Figure 5 —— 51/build_figures.py Figure 5 段（含 panel d 原始基因箭头）
# --------------------------------------------------------------------------------------
def fig5_dryland_context(root, out):
    plt = _rc()
    U = root / "data/upstream"
    D = root / "data/derived"
    cx = _tsv(D / "system_context_summary.tsv")
    fig = plt.figure(figsize=(10.5, 7), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1], width_ratios=[1, 1.18])
    a = fig.add_subplot(gs[0, 0]); c = fig.add_subplot(gs[0, 1])
    e = fig.add_subplot(gs[1, 0]); g = fig.add_subplot(gs[1, 1])
    ss = _tsv(D / "system_state_summary.tsv")
    ss = ss[ss.cohort == "new_desert_strict"].set_index("system").loc[[s for s in SYSTEMS if s != "EctD"]]
    x = np.arange(5)
    a.barh(x, ss.panel_complete, color=TEAL, label="All specified markers")
    a.barh(x, ss.partial, left=ss.panel_complete, color="#D4D8D9", label="Partial marker set")
    for i, row in enumerate(ss.itertuples()):
        a.text(row.panel_complete + row.partial + 1, i, f"{row.panel_complete}/{row.panel_complete + row.partial}",
               va="center", fontsize=9)
    a.set(yticks=x, yticklabels=ss.index, xlabel="Genomes with any marker; labels = all / any", xlim=(0, 95))
    a.invert_yaxis(); _title(a, "a", "Strict dryland marker sets (n = 154)")
    a.legend(frameon=False, fontsize=9, loc="center right", bbox_to_anchor=(1, .39))
    cc = cx[(cx.cohort == "new_desert_strict") & (cx.metric == "compact_10kb")].set_index("system").loc[ss.index]
    for i, row in enumerate(cc.itertuples()):
        c.plot([row.conditional_q025, row.conditional_q975], [i, i], color=TEAL, lw=2)
        c.plot(row.fraction, i, "o", color=TEAL)
        c.text(1.035, i, f"{row.n}/{row.denominator}", va="center", fontsize=10)
    c.set(yticks=x, yticklabels=ss.index, xlabel="Fraction among co-detected genomes", xlim=(0, 1.26))
    c.invert_yaxis(); _title(c, "b", "Same strand and contig within 10 kb")
    support = _tsv(D / "lineage_group_support.tsv")
    ranks = ["class", "order", "family"]
    y = np.arange(3)
    for j, cohort in enumerate(["new_desert_primary", "new_desert_strict"]):
        s = support[support.cohort == cohort].set_index("rank").loc[ranks]
        e.barh(y + (j - .5) * .3, s.n, height=.27, color=TEAL if j == 0 else ORANGE,
               label="Primary" if j == 0 else "Strict")
        for i, n in enumerate(s.n):
            e.text(n + 8, i + (j - .5) * .3, str(n), va="center", fontsize=10)
    e.set(yticks=y, yticklabels=["Class", "Order", "Family"], xlabel="Representatives in supported lineages",
          xlim=(0, 675))
    e.invert_yaxis(); _title(e, "c", "Cross-site lineage support")
    e.legend(frameon=False, fontsize=9, loc="lower right")
    loci = _tsv(U / "context/strict_compact_locus_candidates.tsv")
    selected = []
    cols = ["#176B80", "#DEA158", "#A67EAB"]
    genes = {"K01546": "kdpA", "K01547": "kdpB", "K01548": "kdpC", "K02000": "proV", "K02001": "proW",
             "K02002": "proX", "K06718": "ectA", "K00836": "ectB", "K06720": "ectC",
             "K00108": "betA", "K00130": "betB"}
    for i, system in enumerate(["KdpABC", "ProVWX", "EctABC", "BetAB"]):
        candidates = loci[loci.system == system]
        gid = sorted(candidates.genome_id.unique())[0]
        d = candidates[candidates.genome_id == gid].sort_values("start")
        selected.append(d)
        offset = d.start.min()
        for j, row in enumerate(d.itertuples()):
            st = (row.start - offset) / 1000
            en = (row.end - offset + 1) / 1000
            sign = 1 if row.strand == "+" else -1
            beg = st if sign == 1 else en
            g.add_patch(FancyArrow(beg, i, sign * (en - st), 0, width=.17, head_width=.26,
                                   head_length=min(.15, (en - st) * .35),
                                   length_includes_head=True, facecolor=cols[j % 3], edgecolor="none"))
            g.text((st + en) / 2, i - .22, genes[row.KO], ha="center", fontsize=9)
        g.text(0, i + .31, gid.replace("desert2025__", "") + " · " + str(d.source_or_site.iloc[0]),
               fontsize=8.5, color=".35")
    g.set(yticks=range(4), yticklabels=["KdpABC", "ProVWX", "EctABC", "BetAB"],
          xlabel="Distance from first marked CDS (kb)", xlim=(-.05, 6.4), ylim=(3.7, -.6))
    _title(g, "d", "Illustrative observed loci")
    g.spines["left"].set_visible(False); g.tick_params(axis="y", length=0)
    return _save(fig, out, "figure5_dryland_context")


# --------------------------------------------------------------------------------------
# Figure S1 / S2 —— 34/build_figures.py（conceptual framework / Bayesian profiles）
# --------------------------------------------------------------------------------------
def figs1_framework(root, out):
    plt = _rc(**{"font.family": "DejaVu Sans", "font.size": 10})
    fig, ax = plt.subplots(figsize=(10, 5.7)); ax.set(xlim=(0, 10), ylim=(0, 5.7)); ax.axis("off")

    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=.10", fc=color, ec="#9AABB5", lw=.8))
        ax.text(x + .13, y + h - .22, title, fontsize=11, weight="bold", va="top", color="#173647")
        ax.text(x + .13, y + h - .62, body, fontsize=9.5, va="top", linespacing=1.45, color="#263942")

    box(.15, 3.5, 4.4, 1.6, "CONSTRAINTS",
        "Hydration • salinity / ion chemistry • pH\nTemperature • exposure / shielding\n"
        "Preserve measured units and sampling scale", "#EAF1F7")
    box(5.35, 3.5, 4.4, 1.6, "RESOURCE CONTEXT",
        "Organic-carbon stocks and composition\nEnergy donors / acceptors and mineral context\n"
        "Stocks and elements do not establish usable flux", "#F6F0DE")
    box(.15, .98, 9.6, 1.4, "GENOMIC STRATEGY COMBINATIONS",
        "Cellular maintenance  |  Energy acquisition  |  Carbon assimilation\n"
        "Lineage and detection as competing explanations; strategies need not be mutually exclusive", "#EAF3EE")
    for x in [2.4, 7.6]:
        ax.annotate("", xy=(x, 2.48), xytext=(x, 3.36),
                    arrowprops={"arrowstyle": "->", "color": "#46616D", "lw": 1.5, "linestyle": "--"})
    ax.text(5, 2.92, "Conditional relationships to test", ha="center", fontsize=10, color="#5B6266",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 2})
    ax.text(5, .42, "Bounded implications: geological targeting for life detection | candidate functions "
                    "for engineered resource use", ha="center", fontsize=9)
    ax.text(.1, 5.5, "Conceptual framework — no causal effect or Mars performance estimated",
            fontsize=11, weight="bold")
    return _save(fig, out, "figureS1_framework", dpi=220)


def figs2_bayesian_profiles(root, out):
    plt = _rc(**{"font.family": "DejaVu Sans", "font.size": 10})
    vals = [("Main", 17.70, 14.12, 21.16), ("Shared families", 17.04, 12.66, 21.22),
            ("Prior half", 17.14, 13.61, 20.65), ("Prior double", 17.93, 14.33, 21.42),
            ("Diffuse baseline", 17.43, 14, 20.77)]
    fig, ax = plt.subplots(figsize=(8.5, 4.2)); fig.subplots_adjust(left=.25, right=.98, bottom=.2, top=.86)
    for i, (lab, m, l, u) in enumerate(vals):
        ax.errorbar(m, i, xerr=[[m - l], [u - m]], fmt="o", color="#286B84", capsize=3)
        ax.text(22.1, i, f"{m:.2f} [{l:.2f}, {u:.2f}]", va="center", fontsize=9)
    ax.axvline(0, c="0.65", lw=1)
    ax.set(yticks=range(5), yticklabels=[v[0] for v in vals], xlim=(-1, 33),
           xlabel="Stordalen − Qaidam osmotic-marker support (percentage points)")
    ax.invert_yaxis()
    ax.set_title("Existing Bayesian source contrast · median and 95% credible interval",
                 loc="left", weight="bold", fontsize=11)
    return _save(fig, out, "figureS2_bayesian_profiles", dpi=220)


# --------------------------------------------------------------------------------------
# Figure S3 —— 51/build_figures.py KO-null 段（仅 within-source 零模型；含聚类计数面板）
# --------------------------------------------------------------------------------------
def figs3_null_details(root, out):
    plt = _rc()
    U = root / "data/upstream"
    r = pd.read_csv(U / "panel_repair/recurrence/KO_recurrence_null_results.csv")
    r = r[(r.stratum == "primary_ANI95_representatives") & r.variable_trait]
    rn = r[r.null_model == "within_source_prevalence_preserving"]
    fig = plt.figure(figsize=(10, 5.6), layout="constrained")
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.12])
    for j, dom in enumerate(["Bacteria", "Archaea"]):
        ax = fig.add_subplot(gs[0, j])
        d = rn[rn.domain == dom]
        lim = 0
        for axis, col, mark, lbl in [("cell_maintenance", TEAL, "o", "Maintenance"),
                                     ("energy_acquisition", ORANGE, "^", "Energy"),
                                     ("surface_retention", GREY, "s", "Surface")]:
            g = d[d.strategy_axis == axis]
            xx = g.null_mean / g.n_tips
            yy = g.observed_score / g.n_tips
            ax.scatter(xx, yy, s=32, color=col, marker=mark, alpha=.75, label=lbl)
            lim = max(lim, xx.max(), yy.max())
        ax.plot([0, lim * 1.06], [0, lim * 1.06], "--", color=".6", lw=1)
        ax.set(xlabel="Null mean transitions / tip", ylabel="Observed transitions / tip",
               xlim=(-.006, lim * 1.06), ylim=(-.006, lim * 1.06))
        _title(ax, chr(97 + j), dom)
        if j == 0:
            ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax = fig.add_subplot(gs[0, 2])
    summary = []
    for dom in ["Bacteria", "Archaea"]:
        for axis in ["cell_maintenance", "energy_acquisition"]:
            for null in ["global_prevalence_preserving", "within_source_prevalence_preserving"]:
                d = r[(r.domain == dom) & (r.strategy_axis == axis) & (r.null_model == null)]
                summary.append(dict(domain=dom, axis=axis, null=null, n=len(d),
                                    clustered=int((d.p_less_bh < .05).sum())))
    ss = pd.DataFrame(summary)
    for i, (dom, axis) in enumerate([(d, a) for d in ["Bacteria", "Archaea"]
                                     for a in ["cell_maintenance", "energy_acquisition"]]):
        for j, null in enumerate(["global_prevalence_preserving", "within_source_prevalence_preserving"]):
            row = ss[(ss.domain == dom) & (ss.axis == axis) & (ss["null"] == null)].iloc[0]
            yv = i + (j - .5) * .3
            ax.barh(yv, row.clustered / row.n, height=.26, color=TEAL if j else "#BBCDD2")
            ax.text(.02, yv, f"{row.clustered}/{row.n}", va="center", fontsize=9,
                    color="white" if j else "black")
    ax.set(yticks=range(4),
           yticklabels=["Bacteria\nMaintenance", "Bacteria\nEnergy", "Archaea\nMaintenance", "Archaea\nEnergy"],
           xlabel="Fraction clustered (BH q < 0.05)", xlim=(0, 1.02))
    ax.invert_yaxis(); _title(ax, "c", "Clustering persists within source")
    ax.legend(handles=[Patch(color="#BBCDD2", label="Global null"), Patch(color=TEAL, label="Within-source null")],
              frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(.5, -.16))
    return _save(fig, out, "figureS3_null_details")


# --------------------------------------------------------------------------------------
# Figure S4 —— legacy 435 R 脚本的面板布局/配色/标签（matplotlib 等价移植）
#   (A 策略检出 | B 跨库计数) / (C 比对门槛敏感)，tag A/B/C，单研究四库
# --------------------------------------------------------------------------------------
def figs4_independent_occurrence(root, out):
    plt = _rc(**{"font.family": "DejaVu Sans", "font.size": 8})
    plt.rcParams.update({"axes.linewidth": .35, "xtick.major.width": .3, "ytick.major.width": .3})
    d = root / "data/upstream/figs4/core65"
    a = pd.read_csv(d / "Figure5A_strategy_detection.csv")
    b = pd.read_csv(d / "Figure5B_detection_frequency.csv")
    c = pd.read_csv(d / "Figure5C_threshold_breadth.csv")
    assert int(a.represented_KOs.sum()) == 65, "FigS4 core65 represented_KOs 之和应为 65"
    assert int(b.KO_count.sum()) == int(a.represented_KOs.sum())
    assert int(b.loc[b.metagenomes_with_dual_mate_detection > 0, "KO_count"].sum()) == int(a.detected_KOs.sum())
    assert len(c) == 12
    cols = {"cell_maintenance": "#A56779", "energy_acquisition": "#497778", "surface_retention": "#809DA5"}
    ylabels = {"cell_maintenance": "Maintenance", "energy_acquisition": "Energy",
               "surface_retention": "Surface retention"}
    fig = plt.figure(figsize=(183 / 25.4, 148 / 25.4), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.1])
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, :])
    # A: y = rev(names(cols)) == surface_retention, energy_acquisition, cell_maintenance
    order = ["surface_retention", "energy_acquisition", "cell_maintenance"]
    aa = a.set_index("strategy_axis").loc[order]
    axA.barh(np.arange(3), aa.detection_fraction, color=[cols[k] for k in order], height=.5)
    for i, k in enumerate(order):
        axA.text(aa.loc[k, "detection_fraction"] + .018, i,
                 f"{int(aa.loc[k, 'detected_KOs'])}/{int(aa.loc[k, 'represented_KOs'])}",
                 va="center", fontsize=2.8 * 2.2)
    axA.set(yticks=np.arange(3), yticklabels=[ylabels[k] for k in order],
            xlim=(0, 1.3), xticks=[0, .5, 1], xlabel="Fraction detected")
    axA.set_title("Detection by strategy", loc="left", fontsize=9, fontweight="bold")
    # B
    bb = b.sort_values("metagenomes_with_dual_mate_detection")
    axB.bar(bb.metagenomes_with_dual_mate_detection, bb.KO_count, color="#497778", width=.64)
    for xx, yy in zip(bb.metagenomes_with_dual_mate_detection, bb.KO_count):
        axB.text(xx, yy, f"{int(yy)}", ha="center", va="bottom", fontsize=2.8 * 2.2)
    axB.set(xticks=range(5), xlabel="Metagenomes with detection", ylabel="Number of KOs")
    axB.set_ylim(0, bb.KO_count.max() * 1.13)
    axB.set_title("Across-library occurrence", loc="left", fontsize=9, fontweight="bold")
    # C
    samplecols = {"SRR18183467": "#6C628A", "SRR18183468": "#71948B",
                  "SRR18183469": "#C48B70", "SRR18183470": "#A56779"}
    samplelab = {"SRR18183467": "Last Chance Cave", "SRR18183468": "Indian Tunnel - powder",
                 "SRR18183469": "Indian Tunnel - crystals", "SRR18183470": "Hidden Cave"}
    markers = {"SRR18183467": "s", "SRR18183468": "o", "SRR18183469": "^", "SRR18183470": "D"}
    gates = ["broad", "primary", "stringent"]
    c = c.copy()
    for sid, g in c.groupby("sample_id"):
        sid = str(sid)
        g = g.set_index("gate").loc[gates]
        axC.plot(range(3), g.detected_KOs, marker=markers[sid], lw=.55, ms=4.5,
                 color=samplecols[sid], label=samplelab[sid])
    axC.set(xticks=range(3), xticklabels=["Broad", "Primary", "Stringent"],
            xlabel="Alignment gate", ylabel="Detected KOs")
    axC.set_ylim(c.detected_KOs.min() - 1, c.detected_KOs.max() + 1)
    axC.set_yticks(range(48, 65, 2))
    axC.set_title("Alignment-stringency sensitivity", loc="left", fontsize=9, fontweight="bold")
    axC.legend(frameon=False, fontsize=7, loc="center left", bbox_to_anchor=(1.01, .5))
    for gg, tag in [(axA, "A"), (axB, "B"), (axC, "C")]:
        gg.text(-.02, 1.16, tag, transform=gg.transAxes, fontsize=11, fontweight="bold", va="top")
    return _save(fig, out, "figureS4_independent_occurrence")


# --------------------------------------------------------------------------------------
# Figure 6 —— 52/restore_candidate_figure.py（原始矢量 + 已批准标签修正链）
# --------------------------------------------------------------------------------------
def fig6_dpfq008(root, out):
    import os
    import shutil
    import subprocess
    import xml.etree.ElementTree as ET
    src = root / "data/upstream/fig6_dpfq008/Figure3_DPFQ008_lineage_context.svg"
    s = src.read_text(encoding="utf-8")
    changes = {"cytochrome c-like ×2": "cytochrome c-like",
               "Non-random genomic context": "Context in selected carriers",
               "A cross-environment but lineage-structured family": "A lineage-structured candidate across sources",
               "non-Qaidam pairs:": "pairs excluding Q–Q:"}
    for a, b in changes.items():
        assert a in s, f"Fig6 原始矢量缺少待修正标签：{a!r}"
        s = s.replace(a, b)
    root_el = ET.fromstring(s)
    for el in root_el.iter():
        if el.tag.endswith("text") and any(b in (el.text or "") for b in changes.values()):
            el.attrib.pop("textLength", None)
            el.attrib.pop("lengthAdjust", None)
    out.mkdir(parents=True, exist_ok=True)
    svg_out = out / "figure6_DPFQ008.svg"
    ET.register_namespace("", "http://www.w3.org/2000/svg")
    ET.ElementTree(root_el).write(svg_out, encoding="utf-8", xml_declaration=True)
    made = [str(svg_out)]
    exe = shutil.which("rsvg-convert")
    if not exe:
        for p in ("/opt/homebrew/bin/rsvg-convert", "/usr/local/bin/rsvg-convert",
                  "/usr/bin/rsvg-convert", "/snap/bin/rsvg-convert"):
            if os.path.isfile(p) and os.access(p, os.X_OK):
                exe = p
                break
    if exe:
        for args, target in ((["-f", "pdf", "-o", str(out / "figure6_DPFQ008.pdf")], out / "figure6_DPFQ008.pdf"),
                             (["-w", "2600", "-o", str(out / "figure6_DPFQ008.png")], out / "figure6_DPFQ008.png")):
            subprocess.run([exe, *args, str(svg_out)], check=True)
            made.append(str(target))
        print(f"    fig6 label-restore: rsvg-convert={exe}")
    else:
        try:
            import cairosvg
            cairosvg.svg2pdf(url=str(svg_out), write_to=str(out / "figure6_DPFQ008.pdf"))
            cairosvg.svg2png(url=str(svg_out), write_to=str(out / "figure6_DPFQ008.png"), output_width=2600)
            made += [str(out / "figure6_DPFQ008.pdf"), str(out / "figure6_DPFQ008.png")]
            print("    fig6 label-restore: cairosvg")
        except Exception as exc:  # noqa: BLE001
            print(f"    fig6 label-restore: SVG only ({type(exc).__name__})")
    return made


PLOTTERS = {"fig1": fig1_cohort_strategies, "fig2": fig2_lineage_mosaic, "fig3": fig3_portfolio_uncertainty,
            "fig4": fig4_molecular_configurations, "fig5": fig5_dryland_context, "fig6": fig6_dpfq008,
            "figS1": figs1_framework, "figS2": figs2_bayesian_profiles, "figS3": figs3_null_details,
            "figS4": figs4_independent_occurrence}


def figs5_joint_support(root, out):
    import runpy, subprocess, sys
    mod=runpy.run_path(str(root/'scripts/joint/joint_configuration_support.py'))
    mod['run'](root/'data/upstream',out/'joint')
    prefix=out/'figureS5_joint_configuration_support'
    subprocess.run([sys.executable,str(root/'scripts/figures/plot_joint_support.py'),'--input-dir',str(out/'joint'),'--output-prefix',str(prefix)],check=True)
    return [str(prefix)+'.'+ext for ext in ['png','pdf','svg']]

PLOTTERS['figS5']=figs5_joint_support
