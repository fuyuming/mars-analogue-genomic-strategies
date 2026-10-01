from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,pandas as pd,json
H=Path(__file__).resolve().parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host'
M=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox'];COLS=M+['sulfur_chemolithotrophy','trace_gas_energy']
rows=[];validation=[]
def getx(meta):
 z=pd.get_dummies(meta.lineage_id,dtype=float).to_numpy();s=pd.get_dummies(meta.dataset_id,dtype=float).to_numpy()[:,1:];q=meta[['checkm2_completeness','checkm2_contamination']].to_numpy();q=(q-q.mean(0))/q.std(0)
 xr=np.column_stack([z,q]);return xr,np.column_stack([xr,s])
def estimate(y,xr,xf,w=None):
 if w is None:w=np.ones(len(y))
 w=w/w.mean();mu=np.average(y,axis=0,weights=w);sd=np.sqrt(np.average((y-mu)**2,axis=0,weights=w));z=(y-mu)/sd;z=z*np.sqrt(w[:,None]);a=xr*np.sqrt(w[:,None]);b=xf*np.sqrt(w[:,None]);er=z-a@np.linalg.lstsq(a,z,rcond=None)[0];ef=z-b@np.linalg.lstsq(b,z,rcond=None)[0];sr=(er**2).sum(0);sf=(ef**2).sum(0)
 ma=1-sf[:5].sum()/sr[:5].sum();en=1-sf[5:].sum()/sr[5:].sum();return ma,en,ma-en
for branch in ['core65','restored74']:
 d=pd.read_csv(A/f'{branch}_MAG_module_breadth.tsv.gz',sep='\t');pan=pd.read_csv(A/f'{branch}_panel.tsv',sep='\t');ko=pd.read_csv(A/f'{branch}_presence.csv.gz',usecols=['genome_id','KO','present'])
 for strict in [False,True]:
  label=branch+('_strict' if strict else '_primary');meta=pd.read_csv(O/f'{label}_cohort.tsv',sep='\t').set_index('genome_id');y=d.pivot(index='genome_id',columns='module',values='module_breadth').loc[meta.index,COLS].to_numpy();xr,xf=getx(meta)
  for block in sorted(meta.lineage_id.unique()):
   keep=meta.lineage_id!=block;t=meta[keep];a,b=getx(t);ma,en,de=estimate(y[keep],a,b);removed=meta[~keep]
   rows.append({'branch':branch,'strict':strict,'diagnostic':'drop_lineage','removed':block,'n_removed':len(removed),'label':';'.join(removed.collapsed_rank.unique()),'maintenance':ma,'energy':en,'delta':de})
  p=pan[pan.module=='osmotic_desiccation_salt'];ky=ko[ko.KO.isin(p.KO)].pivot(index='genome_id',columns='KO',values='present').loc[meta.index].astype(float)
  assert np.max(np.abs(ky.mean(axis=1).to_numpy()-y[:,3]))<1e-10
  sub={g:ky[list(v.KO)].mean(axis=1).to_numpy() for g,v in p.groupby('submodule')}
  candidates=[('equal_submodule_weight','none',np.column_stack(list(sub.values())).mean(axis=1))]
  candidates += [('drop_osmotic_submodule',g,ky.drop(columns=list(p[p.submodule==g].KO)).mean(axis=1).to_numpy()) for g in sub]
  candidates += [('drop_osmotic_KO',k,ky.drop(columns=k).mean(axis=1).to_numpy()) for k in ky.columns]
  for typ,label0,v in candidates:
   yy=y.copy();yy[:,3]=v;ma,en,de=estimate(yy,xr,xf);rows.append({'branch':branch,'strict':strict,'diagnostic':typ,'removed':label0,'n_removed':0,'label':'exploratory panel sensitivity','maintenance':ma,'energy':en,'delta':de})
  # Validate saved grouped cross-products against a separate direct row-based least-squares fit.
  li=sorted(meta.lineage_id.unique());si=sorted(meta.dataset_id.unique());rng=np.random.default_rng(20261001);lw=rng.exponential(size=(2000,len(li)));sw=rng.exponential(size=(2000,len(si)));lc=np.array([li.index(l) for l in meta.lineage_id]);sc=np.array([si.index(s) for s in meta.dataset_id])
  for scheme,ww in [('lineage_BB',lw[:,lc]),('lineage_source_product_sensitivity',lw[:,lc]*sw[:,sc])]:
   saved=np.load(O/f'{label}_{scheme}_draws.npz');sr=saved['sr'];sf=saved['sf'];deltas=sf[:,5:].sum(1)/sr[:,5:].sum(1)-sf[:,:5].sum(1)/sr[:,:5].sum(1)
   for idx in [0,117,999,1999]:
    observed=estimate(y,xr,xf,ww[idx])[2];err=abs(observed-deltas[idx]);assert err<1e-9
    validation.append({'branch':branch,'strict':strict,'scheme':scheme,'draw':idx,'absolute_delta_error':err})
 print(branch,'done',flush=True)
pd.DataFrame(rows).to_csv(O/'influence_diagnostics.tsv',sep='\t',index=False);pd.DataFrame(validation).to_csv(O/'direct_row_validation.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).groupby(['branch','strict','diagnostic']).delta.agg(['min','max']).to_string())
