# Expanded primary bacterial comparison

This experiment reproduces the completed bacterial-tree grouping, fixed-panel source contrasts, module decomposition and sensitivity analyses from frozen derived inputs. It does not rerun genome assembly, gene prediction, quality assessment, annotation, ANI or tree inference. The input tree contains 4,088 bacteria and completed with UFBoot convergence correlation 0.993 after 700 iterations. Archaeal tree inference is separate and is not reproduced here.

The primary metadata/module matrix includes 4,465 domain-labelled representatives; scripts restrict their models to bacteria. There are 67 unrooted groups and 4,026 bacteria in the 36 groups represented in at least two of 12 catalogues. Catalogues are composite sources, not independent environmental replicates. The identifier alaska_permafrost_reference denotes a mixed legacy soil-reference catalogue and must not be interpreted as exclusively Alaskan. Gene-prediction histories are retained.

Run with the requirements installed, from any working directory:

```sh
python experiments/primary_bacterial_update/prepare_lineage_blocks.py
python experiments/primary_bacterial_update/primary_paired_uncertainty.py
python experiments/primary_bacterial_update/primary_sensitivity.py
python experiments/primary_bacterial_update/primary_paired_uncertainty.py --block-kind order
python experiments/primary_bacterial_update/verify_statistics_direct.py
python experiments/primary_bacterial_update/module_specificity.py
```

Outputs are written under host/. Run in a fresh checkout to preserve existing derived files. Seed is 20261002; 2,000 exponential lineage weights are used. Endpoints share weights within a stratum, not across primary/strict strata with different group definitions. The latter overlap and are not independent validation.

The primary maintenance-minus-energy source partial-R² difference is 0.030106 (95% conditional lineage-weight interval −0.001863 to 0.053474) for core65 and 0.025100 (−0.006416 to 0.049958) for restored74. Both intervals include zero. The osmotic/ion module exceeds all six other fixed modules in pairwise reweighting intervals in both panels; this exploratory decomposition is not a chemical environmental effect, simultaneous-confidence statement, activity or physiological-investment estimate. Full source-specific results and checks are written as TSV files. Source deletion is refitting, not predictive holdout.

The methods retain known technical and sampling heterogeneity. A completed tree does not validate functional annotations. Manuscripts, private server paths, remote logs and raw protein/genome sequences are not part of this experiment.
