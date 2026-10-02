# Completed expanded strict-cohort analysis

This experiment reproduces tree-derived grouping, fixed-panel source contrasts, paired Bayesian-bootstrap and sensitivity analyses, and three figures from **frozen derived inputs**. It does not rerun annotation, raw genome QC, ANI, alignment or IQ-TREE. See ../final_phylogeny for the upstream tree workflow. The primary-quality tree analysis is not completed by this experiment.

Inputs comprise 1,243 strict-quality ANI representatives (1,193 bacteria, 50 archaea), two actual sequence trees, taxonomy/quality/catalogue metadata without server paths, and two frozen module-breadth panels. The analysis retains 1,154 bacteria in 23 groups represented in at least two catalogues. These catalogues are composite sources, not independent environmental replicates. Gene prediction histories differ, as recorded in the metadata. The label `alaska_permafrost_reference` is a legacy internal identifier; the displayed label is **Legacy soil reference**, not an assertion of Alaskan origin.

Run from any working directory using an environment with requirements.txt installed:

```sh
python experiments/strict_tree_update/prepare_lineage_blocks.py
python experiments/strict_tree_update/strict_paired_uncertainty.py
python experiments/strict_tree_update/strict_sensitivity.py
python experiments/strict_tree_update/strict_paired_uncertainty.py --block-kind order
python experiments/strict_tree_update/verify_statistics_direct.py
python experiments/strict_tree_update/plot_completed_strict.py
```

Scripts use their own directory. Derived TSV/NPY files go to host/; PDF/SVG/PNG figures go to figures/. Run in a fresh checkout if preserving an existing output directory. Reproduction uses seed 20261002 and 2,000 weights. Tiny floating point or figure metadata differences can depend on package versions; published verification records the tested versions and input checksums.

The main core65 source partial R² difference is 0.049266 (maintenance minus energy; conditional lineage-weight 95% interval 0.020200–0.089597); restored74 gives 0.043576 (0.014982–0.084959). Removing the solute/ion module makes both intervals include zero. These are conditional reweighting results, not environmental causal effects, independent replication, or a posterior probability of adaptation. Catalogue-deletion checks refit models and are **not predictive holdout validation**. Order-level grouping and panel/weighting sensitivities are exploratory additions after the first contrast.

Grouping is defined by unrooted edge-separable taxa, avoiding arbitrary root-dependent partitioning. The displayed scaffold keeps one representative leaf per group and preserves its path lengths. Heatmaps summarize all group members. Full-tree support and descriptive terminal-branch flags are separately audited; no post-result outlier exclusion is performed. Six archaeal genomes have light-marker hits; all bacterial light-marker values are zero.

See FIGURE_LEGENDS.md for definitions and display limitations. No wet-lab, competition, resource uptake, catalytic specificity, or Mars activity is established by this experiment.
