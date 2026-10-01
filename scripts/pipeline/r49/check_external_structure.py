from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host';N=R/'48_MAG_portfolio_inference_20261001/host/new_catalogue'
outputs=[];comps=[]
for strict in [False,True]:
 st='strict' if strict else 'primary';m=pd.read_csv(O/f'new_desert_{st}_system_states.tsv',sep='\t').set_index('genome_id')
 for sitepolicy in ['representative_site','single_site_ANI_only']:
  mm=m.copy() if sitepolicy=='representative_site' else m[m.ANI_n_sites==1].copy()
  counts=mm.groupby(['class','site_labels']).size().unstack(fill_value=0);keep=counts.index[(counts>=5).sum(axis=1)>=2];mm=mm[mm['class'].isin(keep)]
  if len(mm)<20:continue
  q=mm[['checkm2_completeness','checkm2_contamination']].to_numpy();q=(q-q.mean(0))/q.std(0);xr=np.column_stack([pd.get_dummies(mm['class'],dtype=float),q]);xf=np.column_stack([xr,pd.get_dummies(mm.site_labels,dtype=float).iloc[:,1:]]);r0=np.linalg.matrix_rank(xr);r1=np.linalg.matrix_rank(xf)
  for branch in ['core65','restored74']:
   pan=pd.read_csv(A/f'{branch}_panel.tsv',sep='\t');pan=pan[pan.module=='osmotic_desiccation_salt'];y=pd.read_csv(N/f'{branch}_presence.tsv.gz',sep='\t',index_col=0).loc[mm.index];keys=sorted(pan.submodule.unique());g=np.column_stack([y[list(pan[pan.submodule==k].KO)].mean(axis=1) for k in keys]);counts0=np.array([(pan.submodule==k).sum() for k in keys]);er=g-xr@np.linalg.lstsq(xr,g,rcond=None)[0];ef=g-xf@np.linalg.lstsq(xf,g,rcond=None)[0];cr=er.T@er;cs=cr-ef.T@ef
   for scheme,w in [('KO_equal',counts0/counts0.sum()),('submodule_equal',np.ones(5)/5)]:
    red=w@cr@w;src=w@cs@w;r=src/red
    outputs.append({'strict':strict,'site_assignment':sitepolicy,'branch':branch,'scheme':scheme,'n':len(mm),'classes':mm['class'].nunique(),'sites':mm.site_labels.nunique(),'site_partial_r2':r,'df_adjusted_site_partial_r2':1-(1-r)*(len(mm)-r0)/(len(mm)-r1),'source_diagonal_sum':np.diag(cs)@(w*w),'source_cross_sum':src-np.diag(cs)@(w*w),'source_sse':src,'reduced_sse':red})
    for i in range(5):
     for j in range(i,5):
      f=1 if i==j else 2;comps.append({'strict':strict,'site_assignment':sitepolicy,'branch':branch,'scheme':scheme,'group1':keys[i],'group2':keys[j],'site_component':f*w[i]*w[j]*cs[i,j],'residual_component':f*w[i]*w[j]*cr[i,j]})
pd.DataFrame(outputs).to_csv(O/'external_weighting_structure.tsv',sep='\t',index=False);pd.DataFrame(comps).to_csv(O/'external_covariance_components.tsv',sep='\t',index=False);print(pd.DataFrame(outputs).to_string(index=False))
