from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,json
H=Path(__file__).resolve().parent;O=H/'host';m=pd.read_csv(O/'primary_paired_analysis_cohort.tsv',sep='\t').set_index('tip_id');b=pd.read_csv(O/'primary4465_module_breadth.tsv.gz',sep='\t');cols=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox','sulfur_chemolithotrophy','trace_gas_energy'];Z=pd.get_dummies(m.lineage_id,dtype=float).to_numpy();S=pd.get_dummies(m.dataset_id,dtype=float).to_numpy()[:,1:];Q=m[['completeness','contamination']].to_numpy();Q=(Q-Q.mean(0))/Q.std(0,ddof=1);a=np.column_stack([Z,Q]);x=np.column_stack([a,S]);levels=sorted(m.lineage_id.unique());g=np.array([levels.index(v) for v in m.lineage_id]);draws=np.load(O/'primary_BB_weights.npy');rows=[]
for branch in ['core65','restored74']:
 y=b[b.branch==branch].pivot(index='tip_id',columns='module',values='breadth').loc[m.index,cols].to_numpy();target=np.load(O/f'{branch}_primary_delta_BB.npy')
 for k in [0,7,27,1999]:
  w=draws[k,g];mu=np.average(y,axis=0,weights=w);var=np.average((y-mu)**2,axis=0,weights=w);sy=(y-mu)/np.sqrt(var);root=np.sqrt(w)[:,None];r0=root*sy-root*a@np.linalg.lstsq(root*a,root*sy,rcond=None)[0];r1=root*sy-root*x@np.linalg.lstsq(root*x,root*sy,rcond=None)[0];delta=lambda r,s:(1-(s[:,:5]**2).sum()/(r[:,:5]**2).sum())-(1-(s[:,5:]**2).sum()/(r[:,5:]**2).sum());d=delta(r0,r1);assert abs(d-target[k])<1e-9;rows.append(dict(branch=branch,draw=k,individual_weighted_lstsq_delta=d,grouped_delta=target[k],absolute_error=abs(d-target[k])))
 r0=y-a@np.linalg.lstsq(a,y,rcond=None)[0];r1=y-x@np.linalg.lstsq(x,y,rcond=None)[0];raw=delta(r0,r1);std=delta(r0/y.std(0),r1/y.std(0));rows.append(dict(branch=branch,draw='unscaled_vs_scaled',individual_weighted_lstsq_delta=raw,grouped_delta=std,absolute_error=abs(raw-std)))
pd.DataFrame(rows).to_csv(O/'direct_weighted_regression_checks.tsv',sep='\t',index=False)
light=b[(b.branch=='core65')&(b.module=='light_energy')&(b.breadth>0)].merge(pd.read_csv(O/'primary4465_integrated_metadata.tsv',sep='\t')[['tip_id','dataset_id','expected_domain']],on='tip_id',how='left',validate='one_to_one');light.to_csv(O/'light_marker_detected.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print('light primary:',len(light),'by domain:',light.expected_domain.value_counts().to_dict())
