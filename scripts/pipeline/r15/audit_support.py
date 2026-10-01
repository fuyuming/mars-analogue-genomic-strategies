"""Observed support and candidate-marker coverage; no new association tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

H=Path(__file__).resolve().parent
O=H/'results'/'support_audit';O.mkdir(parents=True,exist_ok=True)
M=json.loads((H/'inputs/metadata.json').read_text())
def save(d,name):d.to_csv(O/name,sep='\t',index=False)
rows=[];design=[]
for cohort in ['primary','strict']:
 d=pd.read_csv(H/f'inputs/qaidam_{cohort}.tsv',sep='\t')
 cfg=M['environment_scales'][cohort]
 for site,g in d.groupby('site_id'):
  for pressure in ['water','EC','TOC']:
   w=g[f'log2_{pressure}_within'];b=g[f'log2_{pressure}_between']
   rows.append(dict(cohort=cohort,site_id=site,pressure=pressure,n_samples=len(g),n_MAGs=g.n_MAGs.sum(),
    within_sd=w.std(ddof=0),within_range=w.max()-w.min(),site_mean_log2=b.iloc[0],
    within_doubling_in_observed_range=bool(w.max()-w.min()>=1)))
 cols=cfg['environment_columns']+cfg['control_columns'];x=d[['z_'+c for c in cols]].to_numpy(float)
 a=np.c_[np.ones(len(x)),x]
 design.append(dict(cohort=cohort,rows=len(a),columns=a.shape[1],rank=np.linalg.matrix_rank(a),condition_number=np.linalg.cond(a)))
 c=pd.DataFrame(x,columns=cols).corr().rename_axis('variable').reset_index();save(c,f'{cohort}_design_correlations.tsv')
save(pd.DataFrame(rows),'pressure_support_by_site.tsv');save(pd.DataFrame(design),'design_rank.tsv')
d=pd.read_csv(H/'inputs/global_strict.tsv',sep='\t')
cross=pd.crosstab(d.family_path,d.dataset_id);cross['n_sources']=(cross>0).sum(axis=1)
save(cross.rename_axis('family_path').reset_index(),'source_family_overlap.tsv')

# These 88 representatives were selected in prior work from four mixed-status
# families; their coverage is not population prevalence or held-out validation.
p=H.parent/'08_configuration_DPFQ008_20260927/results/DPFQ008_analysis_dataset.tsv'
d=pd.read_csv(p,sep='\t');assert len(d)==88 and d.candidate_detected.sum()==38
markers=[c for c in d if c.endswith('_coverage')]
rows=[];counter=[]
for cohort,g in [('four_mixed_families',d),('strict_90_5',d[d.passes_strict_90_5])]:
 for col in markers:
  pos=np.isclose(g[col],1);carrier=g.candidate_detected.astype(bool)
  rows.append(dict(cohort=cohort,marker=col,n=len(g),n_candidate=int(carrier.sum()),
   marker_complete=int(pos.sum()),complete_candidate_detected=int((pos&carrier).sum()),
   complete_candidate_not_detected=int((pos&~carrier).sum()),
   candidate_detected_marker_incomplete=int((~pos&carrier).sum())))
  for _,r in g[pos&~carrier].iterrows():
   counter.append(dict(cohort=cohort,marker=col,genome_id=r.genome_id,family_path=r.family_path,
    genus=r.genus,site_id=r.site_id,completeness=r.checkm2_completeness,contamination=r.checkm2_contamination))
save(pd.DataFrame(rows),'DPFQ008_observed_marker_coverage.tsv')
save(pd.DataFrame(counter),'DPFQ008_marker_counterexample_candidates.tsv')
print(pd.DataFrame(rows).to_string(index=False))
print(pd.DataFrame(design).to_string(index=False))
