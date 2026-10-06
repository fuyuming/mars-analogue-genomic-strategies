from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,json
H=Path(__file__).resolve().parent;O=H/'host'
m=pd.read_csv(O/'primary_paired_analysis_cohort.tsv',sep='\t').set_index('tip_id');b=pd.read_csv(O/'primary4465_module_breadth.tsv.gz',sep='\t')
M=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox','sulfur_chemolithotrophy','trace_gas_energy']
Z=pd.get_dummies(m.lineage_id,dtype=float).to_numpy();S=pd.get_dummies(m.dataset_id,dtype=float).to_numpy()[:,1:];Q=m[['completeness','contamination']].to_numpy();Q=(Q-Q.mean(0))/Q.std(0,ddof=1);A=np.column_stack([Z,Q]);X=np.column_stack([A,S]);rank=A.shape[1]
idx=list(m.groupby('lineage_id',sort=True).indices.values());W=np.load(O/'primary_BB_weights.npy');assert W.shape==(2000,len(idx))
rows=[];contrasts=[];checks=[]
for branch in ['core65','restored74']:
 y=b[b.branch.eq(branch)].pivot(index='tip_id',columns='module',values='breadth').loc[m.index,M].to_numpy();xx=np.array([X[i].T@X[i] for i in idx]);xy=np.array([X[i].T@y[i] for i in idx]);yy=np.array([(y[i]**2).sum(0) for i in idx])
 def fit(w):
  XX=np.einsum('bg,gij->bij',w,xx);XY=np.einsum('bg,gik->bik',w,xy);YY=w@yy
  bf=np.linalg.solve(XX,XY);br=np.linalg.solve(XX[:,:rank,:rank],XY[:,:rank]);sf=YY-np.einsum('bik,bik->bk',bf,XY);sr=YY-np.einsum('bik,bik->bk',br,XY[:,:rank]);assert np.all(sr>0) and np.all(sf>=-1e-6);return 1-sf/sr
 arr=np.concatenate([fit(W[k:k+100]) for k in range(0,2000,100)]);point=fit(np.ones((1,len(idx))))[0]
 # Independent row-weighted least squares check of the grouped sufficient-statistics implementation.
 for k in [0,7,27,1999]:
  w=np.zeros(len(m))
  for j,ix in enumerate(idx):w[ix]=W[k,j]
  rw=np.sqrt(w)[:,None];r0=rw*y-rw*A@np.linalg.lstsq(rw*A,rw*y,rcond=None)[0];r1=rw*y-rw*X@np.linalg.lstsq(rw*X,rw*y,rcond=None)[0];direct=1-(r1*r1).sum(0)/(r0*r0).sum(0);err=np.max(abs(direct-arr[k]));assert err<1e-9;checks.append(dict(branch=branch,draw=k,max_error=err))
 for j,mod in enumerate(M):rows.append(dict(branch=branch,module=mod,point=point[j],median=np.median(arr[:,j]),q025=np.quantile(arr[:,j],.025),q975=np.quantile(arr[:,j],.975)))
 for j,mod in enumerate(M):
  if j==3:continue
  d=arr[:,3]-arr[:,j];contrasts.append(dict(branch=branch,contrast='osmotic_desiccation_salt_minus_'+mod,point=point[3]-point[j],median=np.median(d),q025=np.quantile(d,.025),q975=np.quantile(d,.975),Pr_positive=np.mean(d>0)))
 np.save(O/f'{branch}_module_source_partial_r2_draws.npy',arr)
pd.DataFrame(rows).to_csv(O/'module_source_partial_r2.tsv',sep='\t',index=False);pd.DataFrame(contrasts).to_csv(O/'osmotic_vs_other_modules.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'module_numeric_validation.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(pd.DataFrame(contrasts).to_string(index=False))
