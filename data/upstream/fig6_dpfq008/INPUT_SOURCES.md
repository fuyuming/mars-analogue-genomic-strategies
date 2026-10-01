# FIG6 (DPFQ008) input-source index

Reproducing FIG6 from raw trees/values is possible with the bundled inputs below;
`scripts/pipeline/r52/restore_candidate_figure.py` (and the `plot fig6` CLI) use
the *approved SVG label-restore* chain: the original vector geometry is preserved
and only four text labels are corrected (52-round decisions). The rendered SVG is
byte-identical to the archived round-52 `figure6_DPFQ008.svg`.

## A. label-restore chain (used by `plot fig6`)

| bundled file | role |
|---|---|
| `Figure3_DPFQ008_lineage_context.svg` | original 183×126 mm vector (label geometry unchanged) |
| `figure_contract.md`, `Figure3_legend.md`, `QA_notes.md`, `run_summary.txt`, `sessionInfo.txt` | original figure contract, legend, QA, run summary, R session |
| `Source_Data/*.csv` | the exact source values plotted in panels a/b/c |

## B. original plotting chain (`legacy/dpfq008/256...R`) argument → bundled input

`256.make_ISME_figure3_DPFQ008_lineage_context.R` takes ten arguments. Where the
original input could be located it is bundled; otherwise the status is stated.

| # | 256 argument | bundled input (relative to this dir) | status |
|---|---|---|---|
| 1 | HMM_HITS | `hits/DPFQ008_HMM_hit_architecture.tsv` (63 strict hits; strict/cxxch columns) | bundled |
| 2 | CONTEXT_STRATA | `context/matched_target_control_strata.tsv` | bundled |
| 3 | NULL_FREQUENCY | `context/null_distribution_frequency.tsv` (level MAG/protein, 100 000 permutations) | bundled |
| 4 | NEIGHBORHOODS | `context/expanded_target_neighborhoods.tsv` | bundled |
| 5 | GENE_TREE | `tree/DPFQ008_strict63.treefile` (IQ-TREE2 LG+F+I+G4) | bundled |
| 6 | MAG_MANIFEST | `congruence/DPFQ008_tree_tip_metadata.tsv` (genome_id, gtdb_taxonomy_r226, ani_cluster_95) | bundled (equivalent MAG manifest) |
| 7 | TREE_CONGRUENCE | `congruence/DPFQ008_gene_host_tree_congruence.tsv` (r=0.666, 9999 perms) | bundled |
| 8 | TREE_SENSITIVITY | `Source_Data/Figure3_tree_sensitivity.csv` (non-Qaidam mask r=0.305, 378 pairs) | bundled |
| 9 | CONVENTIONAL_EVIDENCE | `conventional/priority13_conventional_validation.tsv` (integrated_ipr_count, TM segments) | bundled |
| 10 | OUT_DIR | (produced by the script) | — |

Additional reference trees: `congruence/DPFQ008_gene_tree_58_ANI95.treefile`,
`congruence/DPFQ008_host_tree_58_ANI95.treefile`.

## C. Preserved 52-round decisions (not re-derived)

* Corrected labels — `cytochrome c-like ×2` → `cytochrome c-like`;
  `Non-random genomic context` → `Context in selected carriers`;
  `A cross-environment but lineage-structured family` → `A lineage-structured
  candidate across sources`; `non-Qaidam pairs:` → `pairs excluding Q–Q:`.
* No data or geometry change (`no_data_or_geometry_change = true`).

## D. AI probability comparability

The three AI scores for DPFQ008 are retained in
`data/accessions/Supplementary_Table_S3_candidate_evidence.tsv`. They are
**method-specific scores that are not comparable to one another** (ProtNote
0.815, mDeepFRI 0.114, DPFunc rank/empirical value); no cross-method calibration
is claimed. Model checkpoints / InterPro vectors are **not** bundled.
