# Exploratory subfamily refinement of the frozen method pilot

This analysis refines the 1,289 broad-family-hit proteins from the existing 16-assembly pilot (15 ANI95 groups). It does not estimate community competition, resource consumption, or environmental effects. No genomes were selected by ecological outcomes.

Run with Python 3 and pandas:

```sh
python refine.py --output reproduced
```

The included compressed TSV contains all 169,729 domain hits passing the stated gates, including unassigned subfamilies. It reproduces the localization-aware postprocessing, **not** the entire HMM search from genomes. The raw-search hash, model hash and exact mapping source are recorded in `provenance.json`. The complete proteomes and large HMM database are not bundled.

The first diagnostic retained every passing subfamily, giving 9 single-label protein–family pairs among 671 baseline pairs. That deliberately inclusive envelope confounds strongly overlapping alternatives with simultaneous functions. The additional analysis sorts hits by independent-domain E-value, then score, and greedily retains domains that do not overlap an existing better anchor by more than half the shorter alignment. It considers all overlapping alternatives within 0, 10 or 20 bits of each anchor. These are exploratory sensitivity settings, not calibrated posterior probabilities or pre-registered ecological tests. They are not an exact reproduction of the upstream adjacent-only overlap filter.

There are 705 retained anchors on 663 proteins. At the 20-bit tolerance, 62 domains have one mapped label and no unmapped alternatives; 57 remain after excluding CBM and GT. These remain biochemical annotation candidates: host-associated labels, synthesis, cellular remodeling and accessory redox activity cannot automatically be treated as external carbon uptake. Gene context and transport evidence are pending. Most domains remain unassigned or ambiguous. The broad “human milk polysaccharide” label is not retained by this postprocessing; this is not evidence of biochemical inability.

The comparison is confined to hit refinement. It does not measure proteome-wide sensitivity, prove enzyme specificity, establish available resources, or justify CaCo/RPS competition scores. Counts in `expected/` provide reproducibility checks, not biological validation.

Reference inspected: [run_dbcan overlap processing](https://github.com/bcb-unl/run_dbcan/blob/c23f0d08d7e2678feca7b485bed4b3a4ce7b694a/dbcan/process/process_utils.py). Existing HMM provenance is bounded as recorded; no version number was inferred from individual model dates.
