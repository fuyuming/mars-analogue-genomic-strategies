# Resource-niche annotation pilot — exploratory, outside the manuscript results

This finite experiment tests whether a CAZyme-family-to-substrate endpoint can be used reliably in the existing Mars-analogue study. It does not test competition or resource partitioning and does not replace the manuscript's maintenance–energy or AI candidate analyses.

## Reproduce the completed checks

Python 3 with pandas is required. Run from this directory:

```bash
python fetch_upstream.py
python probe_caco.py
python compare_pilot_parsers.py
python audit_mackay_observations.py --coverage inputs/mackay_coverage.tsv --parents inputs/mackay_parents.tsv --quality inputs/quality607.tsv --out host
```

`fetch_upstream.py` downloads four small official source/data files at fixed commits and verifies their exact hashes. It does not import or run CaCo. The two parser scripts extract only the named parser function using Python AST; synthetic fixtures are explicitly labelled non-biological. Upstream downloads are ignored by Git and retain their original authorship/licensing at the linked repository.

Inputs include a completed HMMER domtblout, a frozen 16-genome selection, derived quality estimates for 607 candidate MAGs, and published Mackay mapping values plus NCBI parent counts. Private working-directory strings were removed from HMMER comment lines only; all data rows are unchanged. The receipt preserves both original and sanitized hashes. Reference outputs were generated from the original raw file; a reproduced summary differs only in `raw_sha256`. No complete genomes, proteins, reference databases or manuscripts are included.

The 16 MAGs were selected before seeing CAZyme results: eight from each cohort, ordered by SHA256 of `NEE58-pilot-v1:` plus the genome name, from the completeness≥90%, contamination<5% stratum. These are method-test inputs, not independent ecological replicates or a holdout cohort. The 607 input candidates yielded 556 genomes at completeness≥50%, contamination<10%, including 129 in the strict stratum; ANI/source deduplication is still required.

## Completed annotation and limits

All 36,521 previously predicted proteins were searched together against the official prebuilt dbCAN V11 (699 profiles; SHA256 `755b244c18f8627a7b446b6999eef7f5e02144ab00715f41417d3971820777d2`). V11 was chosen to match the CaCo implementation, not as a claim to use the latest release. The search used HMMER `hmmsearch --cpu 8 -E 1e-5 --domE 1e-5`; receipt arguments document the search. The wider initial threshold retains raw rows for stricter parsing. The combined target database changes E-values relative to individual-genome searches, so this is a same-raw-output parser comparison, not an end-to-end reproduction of the published CaCo pipeline.

The current and archived parsers retain 539 protein–family pairs. Changing only the HMM-coverage denominator gives 653; retaining all passing families gives 740 under sequence gates or 671 under independent-domain gates. Candidate substrate sets change in 10/16 MAGs between current and all-sequence-gate parsing. These are sensitivity results, not proof that every added family is correct. Fifty proteins have overlapping different-family candidates under the independent-domain gate.

Family-to-substrate assignments are also uncertain: 178/280 genome–substrate assignments after excluding CBM-only support depend on a single family. Family labels can encompass different activities; an assigned label is not a measured environmental resource, substrate uptake, activity or competition. No RPS values, pairwise ecological tests or environmental fits are produced here.

The Mackay 451×16 table has 2,205 dashes, exactly 147 single-parent MAGs×15 unmeasured cross-sample combinations. Dashes must not become zero. Only 304 coassembly MAGs have all 16 numeric entries (297 primary, 69 strict); numeric read-mapping values alone do not validate genome breadth or coexistence. No sample/community inference uses the incomplete table.

## Sources

- [CaCo paper](https://doi.org/10.1073/pnas.2526391123) and [official code](https://github.com/celiosantosjr/CaCo), archived commit `db50a80e974b4daedb217cf971a68c6079fc33e7`, inspected current commit `c07dc2a7811cbe899f7271d5fb0661884b196bfa`.
- [Official dbCAN V11 file](https://pro.unl.edu/dbCAN2/download_file.php?file=dbCAN-HMMdb-V11.txt); exact published size and complete HMM records were checked; no official checksum was available to us.
- [HMMER documentation](https://eddylab.org/software/hmmer/CURRENT/Userguide.pdf).
- Ortiz et al., *Multiple energy sources and metabolic strategies sustain microbial diversity in Antarctic desert soils* (PNAS 2021), supplementary dataset S5; [NCBI PRJNA630822](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA630822).
- [Danakil source paper](https://doi.org/10.1038/s41559-024-02505-6), [PRJNA541281](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA541281); catalogue extra/external records remain explicitly marked in the main accession manifests.
- [dbCAN3 method paper](https://doi.org/10.1093/nar/gkad328) discusses family polyspecificity and substrate inference at subfamily/gene-cluster levels.
