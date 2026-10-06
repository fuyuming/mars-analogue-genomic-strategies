"""Publication figures for the v9 bacterial update.

This script is deliberately a plotting-only layer.  All quantitative values are
read from the completed primary/strict host tables; no tree inference or
statistical refitting is performed here.  The primary tree is pruned only for
display to one deterministic representative leaf per one of the 67 fixed
blocks.  The heatmaps and source tracks still use all 4,088 primary bacterial
leaves.
"""

from __future__ import annotations

import copy
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


# Required Nature-figure SVG/PDF text settings: keep text editable.
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.spines.left": False,
        "legend.frameon": False,
        "xtick.color": "#465057",
        "ytick.color": "#465057",
        "axes.labelcolor": "#344047",
        "text.color": "#27343A",
        "savefig.facecolor": "white",
    }
)


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
FIGURES = HERE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

PRIMARY_HOST = HERE / "inputs" / "primary"
STRICT_HOST = HERE / "inputs" / "strict"
STRICT_MODULE = STRICT_HOST / "module_source_partial_r2.tsv"
PRIMARY_TREE = PRIMARY_HOST / "bac120.treefile"


SOURCE_ORDER = [
    "alaska_permafrost_reference",
    "australian_basalt_lava_tubes_bay_2025",
    "mauna_loa_lava_tube_fishman_2023",
    "qaidam_basin_mmag_1773",
    "stordalen_mire_2019_hybrid_mags",
    "atacama_halite_rainfall_prjna484015",
    "atacama_salt_crust_prjna351262",
    "atacama_TLT_PRJNA1104199",
    "atacama_boulder_fields_PRJNA665391",
    "global_desert_PRJNA832417",
    "mackay_PRJNA630822",
    "danakil_accession_catalogue",
]
SOURCE_LABELS = [
    "Legacy soil reference",
    "Australian lava tubes",
    "Mauna Loa lava tube",
    "Qaidam Basin",
    "Stordalen mire",
    "Atacama halite",
    "Atacama salt crust",
    "Atacama TLT",
    "Atacama boulder fields",
    "Global desert catalogue",
    "Mackay Glacier soils",
    "Danakil catalogue",
]
SOURCE_COLORS = [
    "#335C67",
    "#6D8B74",
    "#A8B86B",
    "#C49A6C",
    "#887F9E",
    "#6F98B8",
    "#9A6E7B",
    "#C17C74",
    "#A77B50",
    "#678D9A",
    "#8C9A73",
    "#B18C72",
]

MODULES = [
    "cold_protein_quality",
    "dna_repair_radiation",
    "dormancy_resuscitation",
    "osmotic_desiccation_salt",
    "oxidative_redox",
    "sulfur_chemolithotrophy",
    "trace_gas_energy",
    "light_energy",
    "biofilm_eps_surface",
]
MODULE_SHORT = [
    "Cold /\nprotein",
    "DNA repair /\nradiation",
    "Dormancy /\nresuscitation",
    "Osmotic /\nion",
    "Oxidative /\nredox",
    "Sulfur\nenergy",
    "Trace-gas\nenergy",
    "Light\nenergy",
    "EPS /\nsurface",
]
MODULES_7 = MODULES[:7]
MODULE_LABELS_7 = [
    "Cold / protein quality",
    "DNA repair / radiation",
    "Dormancy / resuscitation",
    "Osmotic / ion homeostasis",
    "Oxidative / redox",
    "Sulfur chemolithotrophy",
    "Trace-gas energy",
]

PRIMARY_COLOR = "#2F7188"
STRICT_COLOR = "#A9653C"
INK = "#27343A"
MUTED = "#65737A"
GRID = "#E4E9EB"


def save_publication_figure(fig: plt.Figure, basename: str, dpi: int = 450) -> None:
    """Save editable vector files and a high-resolution PNG preview."""

    for ext in ("svg", "pdf"):
        fig.savefig(FIGURES / f"{basename}.{ext}", bbox_inches="tight", pad_inches=0.04)
    fig.savefig(FIGURES / f"{basename}.png", dpi=dpi, bbox_inches="tight", pad_inches=0.04)


def _clean_axis(ax: plt.Axes) -> None:
    for spine in ax.spines.values():
        spine.set_visible(False)


def _short_rank(rank: str) -> str:
    return {"phylum": "p__", "class": "c__", "order": "o__"}.get(rank, "")


def _tree_coordinates(tree):
    """Return depth, y coordinates, and rectangular tree line segments."""

    depths = tree.depths()
    y = {}
    terminals = tree.get_terminals()
    for idx, clade in enumerate(terminals):
        y[clade] = float(idx)
    for clade in tree.get_nonterminals(order="postorder"):
        y[clade] = float(np.mean([y[child] for child in clade.clades]))

    segments = []
    for clade in tree.get_nonterminals(order="preorder"):
        child_y = [y[child] for child in clade.clades]
        segments.append([(depths[clade], min(child_y)), (depths[clade], max(child_y))])
        for child in clade.clades:
            segments.append([(depths[clade], y[child]), (depths[child], y[child])])
    return depths, y, segments


def _read_primary_tree_inputs():
    mapping = pd.read_csv(PRIMARY_HOST / "primary_tree_lineage_mapping.tsv", sep="\t")
    meta = pd.read_csv(PRIMARY_HOST / "primary4465_integrated_metadata.tsv", sep="\t")
    tree = Phylo.read(PRIMARY_TREE, "newick")

    bacteria_meta = meta.loc[meta["domain"].eq("Bacteria")].copy()
    tree_tips = {leaf.name for leaf in tree.get_terminals()}
    mapping_tips = set(mapping["tip_id"])
    if len(tree_tips) != 4088 or len(mapping_tips) != 4088:
        raise ValueError(f"Expected 4,088 primary bacterial tree tips, got {len(tree_tips)} / {len(mapping_tips)}")
    if tree_tips != mapping_tips:
        raise ValueError("Primary tree tips and primary mapping tips are not identical")
    if len(bacteria_meta) != 4088:
        raise ValueError(f"Expected 4,088 primary bacterial metadata rows, got {len(bacteria_meta)}")

    lineage_by_tip = mapping.set_index("tip_id")["lineage_id"].to_dict()
    representative_by_lineage = {}
    # The first leaf encountered in the completed tree is the display-only
    # representative.  This deterministic choice does not alter block means.
    for leaf in tree.get_terminals():
        lineage = lineage_by_tip[leaf.name]
        representative_by_lineage.setdefault(lineage, leaf.name)
    if len(representative_by_lineage) != 67:
        raise ValueError(f"Expected 67 fixed blocks, got {len(representative_by_lineage)}")

    display_tree = copy.deepcopy(tree)
    keep = set(representative_by_lineage.values())
    for leaf in list(display_tree.get_terminals()):
        if leaf.name not in keep:
            display_tree.prune(leaf)
    if len(display_tree.get_terminals()) != 67:
        raise ValueError("Display pruning did not produce exactly 67 representative leaves")

    display_order = [lineage_by_tip[leaf.name] for leaf in display_tree.get_terminals()]
    if len(display_order) != len(set(display_order)):
        raise ValueError("Display tree contains duplicate fixed-block representatives")

    return mapping, bacteria_meta, tree, display_tree, lineage_by_tip, representative_by_lineage, display_order


def plot_figure1() -> None:
    (
        mapping,
        bacteria_meta,
        _full_tree,
        display_tree,
        lineage_by_tip,
        representative_by_lineage,
        display_order,
    ) = _read_primary_tree_inputs()

    block = (
        mapping.groupby("lineage_id", sort=False)
        .agg(
            block_rank=("collapsed_rank", "first"),
            block_name=("collapsed_name", "first"),
            n=("tip_id", "size"),
        )
        .reindex(display_order)
    )
    block["representative_tip_id"] = [representative_by_lineage[x] for x in display_order]
    if int(block["n"].sum()) != 4088:
        raise ValueError("Primary block counts do not sum to 4,088")

    breadth = pd.read_csv(PRIMARY_HOST / "primary4465_module_breadth.tsv.gz", sep="\t")
    breadth = breadth.loc[breadth["branch"].eq("core65")].merge(
        mapping[["tip_id", "lineage_id"]], on="tip_id", how="inner", validate="many_to_one"
    )
    module_means = (
        breadth.groupby(["lineage_id", "module"], sort=False)["breadth"]
        .mean()
        .unstack("module")
        .reindex(index=display_order, columns=MODULES)
    )
    if module_means.isna().any().any() or module_means.shape != (67, 9):
        raise ValueError("Primary core65 module means are incomplete")

    source = mapping[["tip_id", "lineage_id"]].merge(
        bacteria_meta[["tip_id", "dataset_id"]], on="tip_id", how="left", validate="one_to_one"
    )
    if source["dataset_id"].isna().any():
        raise ValueError("Primary source composition is missing dataset IDs")
    source_counts = pd.crosstab(source["lineage_id"], source["dataset_id"]).reindex(
        index=display_order, columns=SOURCE_ORDER, fill_value=0
    )
    source_fractions = source_counts.div(block["n"].to_numpy(), axis=0)
    if int(source_counts.to_numpy().sum()) != 4088:
        raise ValueError("Primary source composition does not sum to 4,088")

    # Figure 1: 67-row display scaffold plus aligned block composition/function tracks.
    n_rows = len(display_order)
    row_y = np.arange(n_rows, dtype=float)
    fig = plt.figure(figsize=(10.2, 15.8))
    gs = fig.add_gridspec(
        1,
        4,
        width_ratios=[4.70, 0.72, 2.65, 3.95],
        left=0.035,
        right=0.985,
        top=0.895,
        bottom=0.245,
        wspace=0.09,
    )
    ax_tree, ax_n, ax_source, ax_heat = [fig.add_subplot(gs[0, i]) for i in range(4)]

    depths, ycoord, segments = _tree_coordinates(display_tree)
    tree_leaf_rows = np.asarray([ycoord[leaf] for leaf in display_tree.get_terminals()], dtype=float)
    expected_rows = np.arange(n_rows, dtype=float)
    if not np.array_equal(tree_leaf_rows, expected_rows):
        raise ValueError("Pruned tree leaf coordinates are not exactly aligned to the 67 display rows")
    ax_tree.add_collection(LineCollection(segments, colors="#64747B", linewidths=0.55, zorder=1))
    max_depth = max(depths[leaf] for leaf in display_tree.get_terminals())
    for leaf in display_tree.get_terminals():
        lineage = lineage_by_tip[leaf.name]
        label_row = display_order.index(lineage)
        label = f"{_short_rank(block.loc[lineage, 'block_rank'])}{block.loc[lineage, 'block_name']}"
        ax_tree.plot(
            [depths[leaf], max_depth * 1.015],
            [ycoord[leaf], ycoord[leaf]],
            color="#C9D1D4",
            lw=0.30,
            ls=":",
            zorder=0,
        )
        ax_tree.text(
            max_depth * 1.035,
            label_row,
            label,
            ha="left",
            va="center",
            fontsize=4.55,
            color=INK,
            clip_on=False,
        )
    ax_tree.set_xlim(-max_depth * 0.02, max_depth * 1.63)
    # Use the exact same y limits as the aligned n/source/heatmap tracks.  The
    # scale bar is drawn in the x-data/axes-fraction transform below so it does
    # not perturb the 67-row alignment.
    ax_tree.set_ylim(n_rows - 0.5, -0.5)
    ax_tree.axis("off")
    scale_y = -0.10
    scale_x = 0.2  # data-space scale bar; label below reports the same value
    ax_tree.plot(
        [0, scale_x],
        [scale_y, scale_y],
        color=INK,
        lw=0.85,
        clip_on=False,
        transform=ax_tree.get_xaxis_transform(),
    )
    ax_tree.text(
        0,
        scale_y - 0.02,
        "0.2 substitutions/site",
        fontsize=5.8,
        va="top",
        clip_on=False,
        transform=ax_tree.get_xaxis_transform(),
    )
    ax_tree.set_title("a  Unrooted bac120 display scaffold", loc="left", fontsize=8.2, fontweight="bold", pad=5)

    ax_n.set_xlim(0, 1)
    ax_n.set_ylim(n_rows - 0.5, -0.5)
    ax_n.axis("off")
    for y, n in zip(row_y, block["n"].astype(int)):
        ax_n.text(0.5, y, f"{n:,}", ha="center", va="center", fontsize=5.2, color=INK)
    ax_n.set_title("b  n", loc="left", fontsize=8.2, fontweight="bold", pad=5)

    ax_source.set_xlim(0, 1)
    ax_source.set_ylim(n_rows - 0.5, -0.5)
    left = np.zeros(n_rows)
    for j, (label, color) in enumerate(zip(SOURCE_LABELS, SOURCE_COLORS)):
        width = source_fractions.iloc[:, j].to_numpy(dtype=float)
        ax_source.barh(
            row_y,
            width,
            left=left,
            height=0.88,
            color=color,
            edgecolor="white",
            linewidth=0.12,
            label=f"S{j + 1}  {label}",
        )
        left += width
    ax_source.set_xticks([0, 0.5, 1.0])
    ax_source.set_xticklabels(["0", "0.5", "1"], fontsize=5.8)
    ax_source.xaxis.tick_bottom()
    ax_source.tick_params(axis="x", length=2, pad=2)
    ax_source.set_yticks([])
    ax_source.grid(axis="x", color=GRID, lw=0.45, zorder=-1)
    ax_source.set_title("c  12-source composition", loc="left", fontsize=8.2, fontweight="bold", pad=5)
    _clean_axis(ax_source)

    heat = ax_heat.imshow(
        module_means.to_numpy(dtype=float),
        aspect="auto",
        interpolation="nearest",
        vmin=0,
        vmax=1,
        cmap="YlGnBu",
    )
    ax_heat.set_ylim(n_rows - 0.5, -0.5)
    ax_heat.set_yticks([])
    ax_heat.set_xticks(np.arange(len(MODULES)))
    ax_heat.set_xticklabels(MODULE_SHORT, fontsize=5.35, rotation=68, ha="right")
    ax_heat.xaxis.tick_bottom()
    ax_heat.tick_params(axis="x", length=0, pad=2)
    ax_heat.set_xticks(np.arange(-0.5, len(MODULES), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax_heat.grid(which="minor", color="white", lw=0.26)
    ax_heat.tick_params(which="minor", length=0)
    ax_heat.set_title("d  Core65 mean marker breadth", loc="left", fontsize=8.2, fontweight="bold", pad=5)
    _clean_axis(ax_heat)
    cbar = fig.colorbar(heat, ax=ax_heat, orientation="vertical", fraction=0.028, pad=0.018, aspect=35)
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.tick_params(labelsize=5.3, length=2)
    cbar.set_label("Mean breadth\n(0–1)", fontsize=5.6, labelpad=3)
    cbar.outline.set_visible(False)

    fig.text(
        0.035,
        0.965,
        "Figure 1 | Primary bacterial sequence tree and fixed-block functional context",
        ha="left",
        va="bottom",
        fontsize=11.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.035,
        0.944,
        "4,088 primary bacteria · 67 fixed blocks · the tree is pruned to one display representative per block",
        ha="left",
        va="bottom",
        fontsize=7.6,
        color=MUTED,
    )
    source_handles = [Patch(facecolor=color, edgecolor="none", label=f"S{j + 1}  {label}") for j, (label, color) in enumerate(zip(SOURCE_LABELS, SOURCE_COLORS))]
    fig.legend(
        handles=source_handles,
        loc="lower center",
        bbox_to_anchor=(0.51, 0.047),
        ncol=4,
        fontsize=5.8,
        handlelength=0.9,
        handleheight=0.7,
        columnspacing=1.05,
        labelspacing=0.55,
        borderaxespad=0,
    )
    fig.text(
        0.035,
        0.143,
        "Display scaffold is unrooted; branch-support labels are hidden. Block tracks use all leaves, and the representative pruning is visual only.",
        fontsize=6.2,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.035,
        0.130,
        "n = bacterial genomes per fixed block; source fractions sum to 1 within each row; module breadth is the core65 mean.",
        fontsize=6.2,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    out = block.reset_index().rename(columns={"lineage_id": "block_id"})
    out["display_order"] = np.arange(1, n_rows + 1)
    for module in MODULES:
        out[f"core65_mean_{module}"] = module_means[module].to_numpy()
    for j, source_key in enumerate(SOURCE_ORDER):
        prefix = f"S{j + 1:02d}"
        out[f"{prefix}_count"] = source_counts[source_key].to_numpy(dtype=int)
        out[f"{prefix}_fraction"] = source_fractions[source_key].to_numpy(dtype=float)
    out.to_csv(FIGURES / "figure1_primary_tree_source.tsv", sep="\t", index=False, float_format="%.12g")
    save_publication_figure(fig, "figure1_primary_tree", dpi=450)
    plt.close(fig)


def _read_pairwise_tables():
    primary = pd.read_csv(PRIMARY_HOST / "primary_paired_uncertainty.tsv", sep="\t")
    strict = pd.read_csv(STRICT_HOST / "strict_paired_uncertainty.tsv", sep="\t")
    primary["cohort"] = "primary"
    strict["cohort"] = "strict"
    primary["cohort_label"] = "Primary bacteria"
    strict["cohort_label"] = "Strict bacteria"
    pair = pd.concat([primary, strict], ignore_index=True)
    pair["cohort_order"] = pair["cohort"].map({"primary": 0, "strict": 1})
    pair["branch_order"] = pair["branch"].map({"core65": 0, "restored74": 1})
    return pair.sort_values(["branch_order", "cohort_order"]).reset_index(drop=True)


def plot_figure2() -> None:
    pair = _read_pairwise_tables()
    expected = {
        ("primary", "core65"): (4026, 36),
        ("primary", "restored74"): (4026, 36),
        ("strict", "core65"): (1154, 23),
        ("strict", "restored74"): (1154, 23),
    }
    for key, (n_expected, blocks_expected) in expected.items():
        row = pair.loc[(pair["cohort"] == key[0]) & (pair["branch"] == key[1])]
        if len(row) != 1 or int(row.iloc[0]["n"]) != n_expected or int(row.iloc[0]["lineage_blocks"]) != blocks_expected:
            raise ValueError(f"Unexpected paired uncertainty row for {key}")

    fig, ax = plt.subplots(figsize=(7.5, 3.6))
    fig.subplots_adjust(left=0.18, right=0.965, top=0.78, bottom=0.235)
    y_base = {"core65": 0.0, "restored74": 1.0}
    y_offset = {"primary": -0.15, "strict": 0.15}
    colors = {"primary": PRIMARY_COLOR, "strict": STRICT_COLOR}
    labels = {"primary": "Primary bacteria (4,026 / 36 blocks)", "strict": "Strict bacteria (1,154 / 23 blocks)"}
    for _, row in pair.iterrows():
        y = y_base[row["branch"]] + y_offset[row["cohort"]]
        lo = 100 * float(row["delta_q025"])
        hi = 100 * float(row["delta_q975"])
        estimate = 100 * float(row["delta_point"])
        color = colors[row["cohort"]]
        ax.plot([lo, hi], [y, y], color=color, lw=2.15, solid_capstyle="round", zorder=2)
        ax.scatter([estimate], [y], s=34, color=color, edgecolor="white", linewidth=0.65, zorder=3)

    ax.axvline(0, color="#6D777C", lw=0.8, ls=(0, (3, 2)), zorder=0)
    ax.set_xlim(-1.5, 10.5)
    ax.set_ylim(1.48, -0.48)
    ax.set_yticks([0, 1])
    ax.set_yticklabels(["core65", "restored74"], fontsize=7.8)
    ax.set_xlabel("Maintenance − energy source partial R² (percentage points)", fontsize=8.1, labelpad=7)
    ax.set_xticks([-1, 0, 2, 4, 6, 8, 10])
    ax.grid(axis="x", color=GRID, lw=0.55)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.tick_params(axis="x", length=3, labelsize=7)
    _clean_axis(ax)
    ax.spines["bottom"].set_visible(True)
    ax.spines["bottom"].set_color("#879298")
    legend = ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color=colors["primary"], lw=2.0, markerfacecolor=colors["primary"], markeredgecolor="white", markersize=5.5, label=labels["primary"]),
            Line2D([0], [0], marker="o", color=colors["strict"], lw=2.0, markerfacecolor=colors["strict"], markeredgecolor="white", markersize=5.5, label=labels["strict"]),
        ],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.20),
        ncol=2,
        fontsize=6.6,
        handlelength=1.3,
        columnspacing=1.4,
    )
    legend.get_frame().set_visible(False)
    fig.text(
        0.18,
        0.935,
        "Figure 2 | Overall maintenance versus energy-source contrast",
        fontsize=10.8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=INK,
    )
    fig.text(
        0.18,
        0.892,
        "Points: reported point estimates; whiskers: paired Bayesian-bootstrap 95% intervals",
        fontsize=7.2,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.18,
        0.055,
        "Primary and strict are overlapping quality layers, not independent replicates; all estimates are conditional on the observed catalogues.",
        fontsize=6.2,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    source = pair[
        [
            "cohort",
            "cohort_label",
            "branch",
            "n",
            "lineage_blocks",
            "catalogues",
            "maintenance_partial_R2",
            "energy_partial_R2",
            "delta_point",
            "delta_BB_mean",
            "delta_q025",
            "delta_q975",
            "Pr_delta_positive",
            "draws",
            "scope",
        ]
    ].copy()
    source.insert(source.columns.get_loc("delta_point") + 1, "delta_point_percentage_points", source["delta_point"] * 100)
    source.insert(source.columns.get_loc("delta_q025") + 1, "delta_q025_percentage_points", source["delta_q025"] * 100)
    source.insert(source.columns.get_loc("delta_q975") + 1, "delta_q975_percentage_points", source["delta_q975"] * 100)
    source.to_csv(FIGURES / "figure2_quality_contrast_source.tsv", sep="\t", index=False, float_format="%.12g")
    save_publication_figure(fig, "figure2_quality_contrast", dpi=450)
    plt.close(fig)


def plot_figure3() -> None:
    primary = pd.read_csv(PRIMARY_HOST / "module_source_partial_r2.tsv", sep="\t")
    strict = pd.read_csv(STRICT_MODULE, sep="\t")
    primary["cohort"] = "primary"
    strict["cohort"] = "strict"
    primary["cohort_label"] = "Primary bacteria"
    strict["cohort_label"] = "Strict bacteria"
    for df in (primary, strict):
        missing = set(MODULES_7) - set(df["module"])
        if missing or set(df["branch"]) != {"core65", "restored74"}:
            raise ValueError(f"Seven-module table is incomplete: {missing}")
    modules = pd.concat([primary, strict], ignore_index=True)
    n_blocks = {
        "primary": {"n": 4026, "lineage_blocks": 36, "catalogues": 12},
        "strict": {"n": 1154, "lineage_blocks": 23, "catalogues": 12},
    }
    for key, values in n_blocks.items():
        for field, value in values.items():
            modules.loc[modules["cohort"].eq(key), field] = value

    # Keep the fixed pre-specified order in both panels.
    y = np.arange(len(MODULES_7), dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.95), sharey=True)
    fig.subplots_adjust(left=0.265, right=0.965, top=0.77, bottom=0.22, wspace=0.10)
    offsets = {"primary": -0.15, "strict": 0.15}
    colors = {"primary": PRIMARY_COLOR, "strict": STRICT_COLOR}
    for panel_idx, branch in enumerate(["core65", "restored74"]):
        ax = axes[panel_idx]
        for cohort in ("primary", "strict"):
            sub = modules.loc[(modules["branch"] == branch) & (modules["cohort"] == cohort)].copy()
            sub["module_order"] = sub["module"].map({m: i for i, m in enumerate(MODULES_7)})
            sub = sub.sort_values("module_order")
            yy = y + offsets[cohort]
            ax.hlines(yy, sub["q025"], sub["q975"], color=colors[cohort], lw=2.0, zorder=2)
            ax.scatter(sub["median"], yy, s=28, color=colors[cohort], edgecolor="white", linewidth=0.55, zorder=3)
        ax.axvline(0, color="#6D777C", lw=0.8, ls=(0, (3, 2)), zorder=0)
        ax.set_xlim(0, 0.27)
        ax.set_xticks([0, 0.05, 0.10, 0.15, 0.20, 0.25])
        ax.set_xticklabels(["0", "0.05", "0.10", "0.15", "0.20", "0.25"], fontsize=6.7)
        ax.grid(axis="x", color=GRID, lw=0.55)
        ax.set_title(branch, loc="left", fontsize=8.5, fontweight="bold", pad=8)
        ax.set_xlabel("Source partial R² (unitless)", fontsize=7.5, labelpad=6)
        ax.tick_params(axis="x", length=3)
        _clean_axis(ax)
        ax.spines["bottom"].set_visible(True)
        ax.spines["bottom"].set_color("#879298")
        if panel_idx == 0:
            ax.set_yticks(y)
            ax.set_yticklabels(MODULE_LABELS_7, fontsize=6.9)
            ax.tick_params(axis="y", length=0, pad=5)
        else:
            ax.tick_params(axis="y", length=0, labelleft=False)
        ax.set_ylim(len(MODULES_7) - 0.5, -0.5)

    axes[0].legend(
        handles=[
            Line2D([0], [0], marker="o", color=PRIMARY_COLOR, lw=2.0, markerfacecolor=PRIMARY_COLOR, markeredgecolor="white", markersize=5.0, label="Primary bacteria (n=4,026; 36 blocks)"),
            Line2D([0], [0], marker="o", color=STRICT_COLOR, lw=2.0, markerfacecolor=STRICT_COLOR, markeredgecolor="white", markersize=5.0, label="Strict bacteria (n=1,154; 23 blocks)"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.27),
        ncol=2,
        fontsize=6.2,
        handlelength=1.2,
        columnspacing=1.25,
    )
    fig.text(
        0.035,
        0.945,
        "Figure 3 | Seven-module source partial R² across bacterial quality layers",
        fontsize=10.8,
        fontweight="bold",
        ha="left",
        va="bottom",
        color=INK,
    )
    fig.text(
        0.035,
        0.907,
        "Dots: reported medians; whiskers: lineage-block Bayesian-bootstrap 95% intervals; all seven modules retained",
        fontsize=7.0,
        color=MUTED,
        ha="left",
        va="bottom",
    )
    fig.text(
        0.035,
        0.072,
        "Partial R² is unitless and conditional on the observed catalogues; the zero line is retained for scale and interpretation.",
        fontsize=6.2,
        color=MUTED,
        ha="left",
        va="bottom",
    )

    # The primary seven-module source table contains point/median/intervals;
    # the strict table additionally carries Pr_positive/scope.  Keep this
    # per-figure source table limited to fields present for both cohorts.
    source = modules[
        [
            "cohort",
            "cohort_label",
            "branch",
            "module",
            "point",
            "median",
            "q025",
            "q975",
            "n",
            "lineage_blocks",
            "catalogues",
        ]
    ].copy()
    source.to_csv(FIGURES / "figure3_module_comparison_source.tsv", sep="\t", index=False, float_format="%.12g")
    save_publication_figure(fig, "figure3_module_comparison", dpi=450)
    plt.close(fig)


def main() -> None:
    plot_figure1()
    plot_figure2()
    plot_figure3()
    print(f"Wrote figures and source tables to {FIGURES}")


if __name__ == "__main__":
    main()
