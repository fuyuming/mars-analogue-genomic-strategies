# Joint-configuration support inputs (portable CLI)

`scripts/joint/joint_configuration_support.py` is the publication-facing,
path-parameterised version of the frozen diagnostic described in
`analysis_contract.md`. It reads **only** the small bundled tables below; no
private filesystem path appears in the script.

| CLI role | bundled path (relative to `data/upstream/`) | role |
|---|---|---|
| original cohort presence | `panel_repair/amended_inputs/core65_presence.csv.gz` | per-genome KO presence + family/genus, 4,213-MAG released cohort |
| new dryland cohort presence | `cohorts/new_catalogue/restored74_presence.tsv.gz` | per-genome KO presence, added dryland catalogue |
| per-cohort system states | `cohorts/{original_primary,original_strict,new_desert_primary,new_desert_strict}_system_states.tsv` | panel_complete / partial / not_detected calls |
| representative metadata | `cohorts/new_desert_{primary,strict}_independent_representatives.tsv` | family + source/site unit labels for the new cohort |

Reference outputs produced by the authoring host are stored under
`reference_outputs/` for comparison (identical schema).

## Frozen rules (do not alter)

* Four pairs: `KdpABC`, `ProVWX`, `EctABC`, `BetAB`, each against co-detection of
  `K03518`/`K03519`/`K03520` (Cox/Cut-family markers). Co-detection does **not**
  establish aerobic/atmospheric CO oxidation — sequence validation is required.
* Joint-state support gate: all four binary states have `>= gate` representative
  MAGs within a family **and** `>= 2` source/site groups each with `>= 3` MAGs.
  Declared thresholds: **2 / 3 / 5**.
* Descriptive only: no p-values, no environment model, no causal/synergy claim.
* "Not all detected" is not biological absence; missing family assignments are
  excluded from family conditioning and counted separately.
