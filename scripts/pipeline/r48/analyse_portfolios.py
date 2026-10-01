"""Paired block Bayesian bootstrap of existing corrected MAG portfolio estimand."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,hashlib,itertools
import numpy as np
import pandas as pd
H=Path(__file__).resolve().parent; ROOT=H.parents[1]; RAW=ROOT.parent/'code/result_raw'
A=RAW/'KO_panel_repair_20260926/amended_inputs'; OUT=H/'host';OUT.mkdir(exist_ok=True)
LP=RAW/'ISME_Figure2_monophyletic_lineage_source_20260717_v2/genome_to_monophyletic_lineage.csv.gz'
SEED=20261001; B=2000
M=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox']
E=['sulfur_chemolithotrophy','trace_gas_energy']; COLS=M+E
manifest=[]; checks=[]; results=[]; module_rows=[]; panels=[]
def read(path,**kw):
 manifest.append({'path':str(path.resolve()),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
 return pd.read_csv(path,**kw)
L=read(LP)
def prepare(branch,strict):
 d=read(A/f'{branch}_MAG_module_breadth.tsv.gz',sep='\t')
 d=d[d.species_representative & (d.gtdb_domain=='Bacteria') & ((not strict) | d.passes_strict_90_5)]
 meta=d[['genome_id','dataset_id','checkm2_completeness','checkm2_contamination']].drop_duplicates().merge(L[['genome_id','lineage_id','collapsed_rank']],on='genome_id',validate='one_to_one')
 keep=meta.groupby('lineage_id').dataset_id.nunique();meta=meta[meta.lineage_id.isin(keep[keep>=2].index)].sort_values('genome_id').set_index('genome_id')
 y=d.pivot(index='genome_id',columns='module',values='module_breadth').loc[meta.index,COLS].to_numpy(float)
 assert not np.isnan(y).any()
 assert d[d.module=='light_energy'].module_breadth.var()==0
 Z=pd.get_dummies(meta.lineage_id,dtype=float).to_numpy(); S=pd.get_dummies(meta.dataset_id,dtype=float).to_numpy()[:,1:]
 quality=meta[['checkm2_completeness','checkm2_contamination']].to_numpy();quality=(quality-quality.mean(0))/quality.std(0,ddof=1)
 xr=np.column_stack([Z,quality]);x=np.column_stack([xr,S]);r=xr.shape[1]
 assert np.linalg.matrix_rank(x)==x.shape[1]
 cells=meta.reset_index().groupby(['lineage_id','dataset_id'],sort=True).indices
 C=list(cells); ngr=len(C);xx=[];xy=[];yy=[];ys=[];ns=[]
 for c in C:
  idx=cells[c];xx.append(x[idx].T@x[idx]);xy.append(x[idx].T@y[idx]);yy.append((y[idx]**2).sum(0));ys.append(y[idx].sum(0));ns.append(len(idx))
 stats=[np.array(z) for z in (xx,xy,yy,ys,ns)]
 meta.to_csv(OUT/f'{branch}_{"strict" if strict else "primary"}_cohort.tsv',sep='\t')
 return meta,y,x,r,C,stats

def fit_weights(w,stats,r):
 xx,xy,yy,ys,ns=stats
 xx=np.einsum('bg,gij->bij',w,xx);xy=np.einsum('bg,gik->bik',w,xy)
 yy=w@yy;ys=w@ys;n=w@ns
 var=yy/n[:,None]-(ys/n[:,None])**2
 assert np.all(var>0)
 bf=np.linalg.solve(xx,xy);br=np.linalg.solve(xx[:,:r,:r],xy[:,:r])
 ssef=yy-np.einsum('bik,bik->bk',bf,xy)
 sser=yy-np.einsum('bik,bik->bk',br,xy[:,:r])
 assert np.all(ssef>=-1e-8) and np.all(sser>=ssef-1e-8)
 return sser/var,ssef/var

def part(sr,sf,ix):return 1-sf[:,ix].sum(1)/sr[:,ix].sum(1)
def summary(branch,strict,scheme,sr,sf,name,mi,ei):
 a=part(sr,sf,mi);b=part(sr,sf,ei);delta=a-b
 row={'branch':branch,'strict':strict,'scheme':scheme,'comparison':name,'draws':len(delta),'maintenance_mean':a.mean(),'energy_mean':b.mean(),'delta_mean':delta.mean(),'delta_median':np.median(delta),'delta_q025':np.quantile(delta,.025),'delta_q975':np.quantile(delta,.975),'prob_delta_positive':np.mean(delta>0)}
 results.append(row)
 return delta
def main():
 for branch in ['core65','restored74']:
  for strict in [False,True]:
   label=branch+('_strict' if strict else '_primary');meta,y,x,r,cells,stats=prepare(branch,strict)
   sr0,sf0=fit_weights(np.ones((1,len(cells))),stats,r)
   # Independent individual-level QR reconstruction validates grouped sufficient stats.
   z=(y-y.mean(0))/y.std(0,ddof=1)
   rr=z-x[:,:r]@np.linalg.lstsq(x[:,:r],z,rcond=None)[0]
   rf=z-x@np.linalg.lstsq(x,z,rcond=None)[0]
   expected=read(RAW/f'KO_panel_repair_20260926/portfolio_{branch}/portfolio_source_partition_summary.tsv',sep='\t')
   ref=expected[expected.analysis==('strict_bacteria_tree_lineage' if strict else 'primary_bacteria_tree_lineage')].iloc[0]
   ma=part(sr0,sf0,list(range(5)))[0];en=part(sr0,sf0,[5,6])[0]
   assert abs(ma-ref.maintenance_source_partial_r2)<1e-10
   assert abs(en-ref.energy_source_partial_r2)<1e-10
   assert abs((1-(rf[:,:5]**2).sum()/(rr[:,:5]**2).sum())-ma)<1e-10
   checks.append({'branch':branch,'strict':strict,'n':len(meta),'blocks':meta.lineage_id.nunique(),'sources':meta.dataset_id.nunique(),'maintenance':ma,'energy':en,'delta':ma-en,'rank_reduced':r,'rank_full':x.shape[1],'max_reference_error':max(abs(ma-ref.maintenance_source_partial_r2),abs(en-ref.energy_source_partial_r2))})
   for j,col in enumerate(COLS):
    module_rows.append({'branch':branch,'strict':strict,'module':col,'source_partial_r2':float(1-sf0[0,j]/sr0[0,j]),'breadth_mean':y[:,j].mean(),'breadth_sd':y[:,j].std(ddof=1),'reduced_standardized_sse':sr0[0,j],'portfolio_residual_weight':sr0[0,j]/sr0[0,:5 if j<5 else 7].sum() if j<5 else sr0[0,j]/sr0[0,5:].sum()})
   # All panel comparisons specified before new results.
   comparisons=[('full',list(range(5)),[5,6])]
   comparisons += [('drop_'+M[j],[k for k in range(5) if k!=j],[5,6]) for j in range(5)]
   comparisons += [('two_maintenance_'+'+'.join(M[j] for j in ix),list(ix),[5,6]) for ix in itertools.combinations(range(5),2)]
   comparisons += [('energy_only_'+E[j],list(range(5)),[5+j]) for j in range(2)]
   for name,mi,ei in comparisons:
    panels.append({'branch':branch,'strict':strict,'comparison':name,'maintenance_partial_r2':part(sr0,sf0,mi)[0],'energy_partial_r2':part(sr0,sf0,ei)[0],'delta':(part(sr0,sf0,mi)-part(sr0,sf0,ei))[0]})
   for ei in itertools.combinations(range(7),2):
    mi=[j for j in range(7) if j not in ei]
    panels.append({'branch':branch,'strict':strict,'comparison':'relabel_energy_'+'+'.join(COLS[j] for j in ei),'maintenance_partial_r2':part(sr0,sf0,mi)[0],'energy_partial_r2':part(sr0,sf0,list(ei))[0],'delta':(part(sr0,sf0,mi)-part(sr0,sf0,list(ei)))[0]})
   li=sorted(meta.lineage_id.unique());si=sorted(meta.dataset_id.unique());lc=np.array([li.index(c[0]) for c in cells]);sc=np.array([si.index(c[1]) for c in cells])
   rng=np.random.default_rng(SEED)
   lw=rng.exponential(size=(B,len(li)));sw=rng.exponential(size=(B,len(si)))
   for scheme,w in [('lineage_BB',lw[:,lc]),('lineage_source_product_sensitivity',lw[:,lc]*sw[:,sc])]:
    chunks=[fit_weights(w[i:i+100],stats,r) for i in range(0,B,100)]
    sr=np.concatenate([c[0] for c in chunks]);sf=np.concatenate([c[1] for c in chunks])
    delta=summary(branch,strict,scheme,sr,sf,'full',list(range(5)),[5,6])
    for name,mi,ei in comparisons[1:]:summary(branch,strict,scheme,sr,sf,name,mi,ei)
    np.savez_compressed(OUT/f'{label}_{scheme}_draws.npz',sr=sr,sf=sf,modules=COLS)
   print(label,'n',len(meta),'delta',ma-en,flush=True)
 pd.DataFrame(checks).to_csv(OUT/'point_reproduction.tsv',sep='\t',index=False)
 pd.DataFrame(results).to_csv(OUT/'paired_uncertainty.tsv',sep='\t',index=False)
 pd.DataFrame(module_rows).to_csv(OUT/'module_decomposition.tsv',sep='\t',index=False)
 pd.DataFrame(panels).to_csv(OUT/'all_panel_comparisons.tsv',sep='\t',index=False)
 (OUT/'input_manifest.json').write_text(json.dumps(list({p['path']:p for p in manifest}.values()),indent=2))
 print(pd.DataFrame(results).query("comparison=='full'").to_string(index=False))

if __name__=='__main__':main()
