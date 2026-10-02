# Resource-context and new-cohort endpoint audit (round 61)

These diagnostics reuse completed uniform Kofam calls; they do not rerun annotations, claim uptake, or fit a new ecological mechanism. Run with Python 3:

```bash
python experiments/resource_context/audit_context.py
python experiments/resource_context/audit_endpoint_support.py
```

Small frozen inputs are included. The context table has 57 targets and 431 distinct window proteins **including the targets**: 374 distinct non-target proteins, 134 with accepted KO labels. Nine targets have alternative KO annotations requiring adjudication before treating subfamily substrate labels as external carbon use. Generic transporter text matches do not establish substrate specificity, direction, or a complete uptake pathway. Two local sugar-enzyme/transport neighbourhoods warrant follow-up but are not validated resource-acquisition systems.

The 556-genome x 74-KO table is extracted from the completed added-cohort Kofam mapper, whose SHA is in export_receipt.json. The 74-KO panel and marker-combination definitions are inherited from earlier audited analyses. The scripts select the frozen quality-specific global ANI representatives using ../final_phylogeny/frozen_genome_manifest.tsv; primary/new Mackay=386 and Danakil=113, strict=107 and16. Same-family four-state gate checks describe support, not hypothesis tests. Resource-unit checks apply to source-verified Mackay records; Danakil sample-chemistry eligibility is explicitly not evaluated by that check.

Detected labels and reference family assignments are uncertain biological annotations. A zero denotes no accepted annotation under the fixed pipeline, not biological absence. Quality filtering and domain imbalance must be examined before ecological interpretation. These outputs do not replace the manuscript figures or establish competitive interactions.
