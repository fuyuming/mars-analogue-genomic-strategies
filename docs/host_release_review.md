# Initial public-code review, 1 October 2026

The host ran the reviewed code in a new directory with only bundled inputs. Input checks, all 11 figure renderings and the joint-configuration diagnostic completed successfully. Four joint TSV outputs match the host numerical results byte for byte; the Fig6 SVG matches the archived corrected figure byte for byte. All 11 outputs were visually compared for panel content, axes and labels. Raster byte identity is not claimed for re-rendered figures.

The first wrapper implementation failed this review: it omitted the Fig2 trees, misread metadata as numeric data in Fig4, conflated labels for two nulls in FigS3 and simplified other panels. It was replaced with adaptations of the original plotting code before public release. Earlier local commits are excluded from the public branch.

This verifies the small-input reproduction layer, not a fresh run from raw reads. Heavy legacy analyses are archived and clearly identified as not rerun here. The added-cohort workflow has been launched in its download-wait stage; uniform quality, annotation, taxonomy, ANI and ecological admission results were not yet available at this release. No manuscript or private execution logs are included.
