# Expanded-cohort phylogeny (round 60)

**Status: marker extraction launched on 2026-10-02; new trees and revised ecological results are not yet available.** This workflow does not replace the currently published manuscript figures until inference and biological review finish.

The frozen manifest contains 4,477 unique genomes: 4,466 primary representatives and 1,243 strict representatives, with 1,232 shared. Quality-specific ANI representatives differ; strict trees are inferred separately, not pruned from the primary tree. Marker filtering may reduce final tip counts, with explicit exclusion records. Source names are catalogue identifiers, not ecological replicates. In particular, the historical `alaska_permafrost_reference` identifier is displayed as **Legacy soil reference**; it does not assign Alaska geography. Environmental-link conflicts stay in the genomic trees but are ineligible for those environmental associations.

Requirements: Python 3, GTDB-Tk 2.5.2 with r226 reference package, and IQ-TREE 2.0.7. External sequence files and databases are not bundled. Earlier analysis inputs, accession tables and source audits remain necessary to reconstruct all non-accession legacy IDs; this is a frozen-input tree workflow, not a complete raw-data reconstruction claim.

1. Edit `config.example.json` to give absolute tool/database paths and a **new** output directory; save as a local config.
2. Prepare a TSV containing exactly `genome_id` and `fasta_path` for every frozen manifest genome. FASTA files must be uncompressed.
3. Bind inputs, then launch:

```bash
python bind_inputs.py --paths /path/to/local_paths.tsv --config /path/to/local_config.json
python run_final_phylogeny.py --config /path/to/fresh/output/config.json
```

The runner records sequence hashes and tool versions, applies the official bac120/ar53 canonical masks and a 10% alignment coverage gate, checks IDs and domains, and infers both quality sets with LG+F+R10, 1,000 ultrafast bootstraps and 1,000 SH-aLRT replicates. Markers use 32 threads. Two quality-set inference jobs run concurrently with 32 threads and a 96 GB memory limit each. These resources and versions are explicit, not automatically inferred from the host.

Completed stages have hash-checked receipts. An attempted stage without a completion receipt **stops for diagnosis**: this script will not blindly restart or overwrite it. Preserve IQ-TREE checkpoints for an explicitly reviewed resume; never add `--redo` to bypass a failure. Tree acceptance also requires normal exit, report/support files and exact tip matching. Branch support, long branches, old/new overlap and manuscript implications require subsequent scientific review.
