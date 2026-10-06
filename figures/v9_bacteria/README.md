# Expanded bacterial manuscript figures

Run `python -m pip install -r requirements.txt`, then `python plot_v9_figures.py` in this directory. Python 3.12 was used. Outputs are PNG/PDF/SVG and source TSV in `figures/`.

This plotting layer uses frozen, public-derived inputs. Statistical reproduction is in `experiments/primary_bacterial_update` and the existing strict/module analysis folders. It does not rerun raw genome processing. Figure 1 uses the actual bacterial sequence tree pruned to one display representative per fixed block; all 4,088 genomes contribute to tracks. No aggregate branch support is inferred. Figures 2–3 describe conditional source structure; quality subsets overlap, source confounds environment/geography/study, and the overall primary contrast includes zero. Module contrasts are exploratory and marginal.

Historical `alaska_permafrost_reference` is displayed as Legacy soil reference because provenance is mixed. No original manuscript or full-text literature is included. Arial is preferred, with free-font fallbacks; font substitution can change visual geometry without changing source data.
