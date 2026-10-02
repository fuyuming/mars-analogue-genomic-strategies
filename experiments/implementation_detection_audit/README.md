# Implementation-detection and comparability audit

This is an exploratory audit of marker scope, local gene organization and sampling support. It does **not** establish alternative physiological strategies, adaptation, competition, environmental causality or activity.

Frozen inputs include the original 74-marker calls (5,477 genomes), final primary/strict ANI representatives, source-checked Mackay environmental links, a separately extracted 10-KO table, and original Prodigal coordinates for selected proteins. Original server paths are removed. Five overlapping KO calls agree with the original panel in all 27,385 comparisons. No annotation or phylogenetic inference is rerun here; original annotation histories remain heterogeneous. The extension does not replace the original manuscript panels.

Run from the repository root with Python, numpy and pandas:

```sh
python experiments/implementation_detection_audit/check_identifiability.py
python experiments/implementation_detection_audit/check_expanded_detection.py
python experiments/implementation_detection_audit/audit_opu_organization.py
python experiments/implementation_detection_audit/compare_transport_signatures.py
python experiments/implementation_detection_audit/check_finer_lineage.py
```

Outputs go to host/. Use a clean checkout to preserve existing outputs. These steps reproduce analysis from frozen derived inputs; they do not claim raw-sequence-to-annotation reproduction. Input hashes are in input_manifest.json.

A marker signature requires all **listed** KOs, not a validated functional pathway. Non-complete signatures combine partial and absent calls and must not be interpreted as true absence. Kdp regulation is separated from pump-subunit signatures. BetAB is a choline-to-betaine branch, not proof of de novo synthesis or the only known route. ProVWX and K05845/6/7 signatures have different annotation scopes and need not have identical substrates. See KEGG M00555 and Transporters ko02000.

For the additional ABC signature, local organization requires three distinct ORFs, same contig and strand, ORF-number span ≤5 and maximum inter-ORF gap ≤2,000 bp. A separate count requires all three original Prodigal partial flags to be 00. Failure of this screen is not proof of functional absence; success is not operon or activity validation. Coordinates are available only for the selected newly added cohorts, not for every catalogue.

Results: 144 added-cohort primary/strict representatives have all three additional ABC KOs; 78 have a local triplet and 66 have a triplet with all ORFs predicted non-partial. Of 100/31 source-checked primary/strict Mackay MAGs, 24/4 have all three KOs; 10/2 retain both local-organization and non-partial criteria. The strict environmental subset is insufficient for the proposed multivariable mechanism model; no such model is fitted.

Coarse screening requires ≥3 genomes per single-signature state, each spanning ≥2 catalogues. At genus level only primary Rubrobacter_F passes; its two groups have median completeness around 65–70%, and no strict genus passes. Consequently single-signature states are not accepted as evidence of functional replacement. Catalogue counts are not independent environmental replicates. The old `alaska_permafrost_reference` identifier denotes a mixed Legacy soil reference catalogue, not a uniform Alaskan environment.
