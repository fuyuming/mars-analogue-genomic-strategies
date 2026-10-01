from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,json,hashlib
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host'
subrows=[];decomp=[];totals=[];coeffs=[]
for branch in ['core65','restored74']:
 pan=pd.read_csv(A/f'{branch}_panel.tsv',sep='\t');pan=pan[pan.module=='osmotic_desiccation_salt'];data=pd.read_csv(A/f'{branch}_presence.csv.gz',usecols=['genome_id','KO','present']);data=data[data.KO.isin(pan.KO)]
 for strict in [False,True]:
  label=branch+('_strict' if strict else '_primary');meta=pd.read_csv(R/f'48_MAG_portfolio_inference_20261001/host/{label}_cohort.tsv',sep='\t').set_index('genome_id')
  ky=data.pivot(index='genome_id',columns='KO',values='present').loc[meta.index].astype(float);keys=sorted(pan.submodule.unique());g=np.column_stack([ky[list(pan[pan.submodule==k].KO)].mean(axis=1) for k in keys]);counts=np.array([(pan.submodule==k).sum() for k in keys]);q=meta[['checkm2_completeness','checkm2_contamination']].to_numpy();q=(q-q.mean(0))/q.std(0);xr=np.column_stack([pd.get_dummies(meta.lineage_id,dtype=float),q]);xf=np.column_stack([xr,pd.get_dummies(meta.dataset_id,dtype=float).iloc[:,1:]])
  er=g-xr@np.linalg.lstsq(xr,g,rcond=None)[0];ef=g-xf@np.linalg.lstsq(xf,g,rcond=None)[0];cr=er.T@er;cf=ef.T@ef;cs=cr-cf
  assert np.linalg.eigvalsh(cs).min()>-1e-9
  for j,k in enumerate(keys):
   subrows.append({'branch':branch,'strict':strict,'submodule':k,'n_KO':counts[j],'mean_breadth':g[:,j].mean(),'var_breadth':g[:,j].var(),'source_partial_r2':cs[j,j]/cr[j,j],'reduced_sse':cr[j,j],'source_sse':cs[j,j]})
  for scheme,w in [('KO_equal',counts/counts.sum()),('submodule_equal',np.ones(5)/5)]:
   red=w@cr@w;src=w@cs@w;frac=src/red
   totals.append({'branch':branch,'strict':strict,'scheme':scheme,'module_source_partial_r2':frac,'source_sse':src,'reduced_sse':red,'source_diagonal_sum':np.diag(cs)@(w*w),'source_cross_sum':src-np.diag(cs)@(w*w),'residual_diagonal_sum':np.diag(cr)@(w*w),'residual_cross_sum':red-np.diag(cr)@(w*w)})
   assert abs(1-((ef@w)**2).sum()/((er@w)**2).sum()-frac)<1e-12
   for i in range(5):
    for j in range(i,5):
     f=1 if i==j else 2
     decomp.append({'branch':branch,'strict':strict,'scheme':scheme,'group1':keys[i],'group2':keys[j],'source_component':f*w[i]*w[j]*cs[i,j],'residual_component':f*w[i]*w[j]*cr[i,j],'source_fraction_component':f*w[i]*w[j]*cs[i,j]/src,'reduced_residual_correlation':cr[i,j]/np.sqrt(cr[i,i]*cr[j,j])})
  # Within-system co-detection correlation is descriptive, conditional on observed genes.
  for k,pp in pan.groupby('submodule'):
   yy=ky[list(pp.KO)].to_numpy();rr=yy-xr@np.linalg.lstsq(xr,yy,rcond=None)[0];co=rr.T@rr;ds=np.sqrt(np.diag(co));cc=co/np.outer(ds,ds)
   for i in range(len(pp)):
    for j in range(i+1,len(pp)):
     coeffs.append({'branch':branch,'strict':strict,'submodule':k,'KO1':pp.KO.iloc[i],'KO2':pp.KO.iloc[j],'quality_lineage_adjusted_correlation':cc[i,j]})
pd.DataFrame(subrows).to_csv(O/'submodule_decomposition.tsv',sep='\t',index=False);pd.DataFrame(totals).to_csv(O/'weighting_matrix_totals.tsv',sep='\t',index=False);pd.DataFrame(decomp).to_csv(O/'covariance_components.tsv',sep='\t',index=False);pd.DataFrame(coeffs).to_csv(O/'within_submodule_residual_correlations.tsv',sep='\t',index=False)
print(pd.DataFrame(subrows).query("branch=='core65' and strict==False").to_string(index=False));print(pd.DataFrame(totals).to_string(index=False))
