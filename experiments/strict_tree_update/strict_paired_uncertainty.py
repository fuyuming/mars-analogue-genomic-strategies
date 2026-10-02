from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import numpy as np,pandas as pd,json,argparse
ap=argparse.ArgumentParser();ap.add_argument('--block-kind',choices=['tree','order'],default='tree');args=ap.parse_args()
H=Path(__file__).resolve().parent;I=H/'host';O=I if args.block_kind=='tree' else I/'order_sensitivity';O.mkdir(exist_ok=True);meta=pd.read_csv(I/'strict1243_integrated_metadata.tsv',sep='\t').set_index('tip_id');lin=pd.read_csv(I/'strict_tree_lineage_mapping.tsv',sep='\t').set_index('tip_id');b=pd.read_csv(I/'strict1243_module_breadth.tsv.gz',sep='\t');M=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox'];E=['sulfur_chemolithotrophy','trace_gas_energy'];cols=M+E;results=[];checks=[];decomp=[]
m=meta.join(lin[['lineage_id']]);m=m[m.expected_domain=='Bacteria']
if args.block_kind=='order':m=m[m['order'].ne('Unclassified')].copy();m['lineage_id']='order:'+m['order']
counts=m.groupby('lineage_id').dataset_id.nunique();m=m[m.lineage_id.isin(counts[counts>=2].index)].sort_index();m.to_csv(O/'strict_paired_analysis_cohort.tsv',sep='\t');Z=pd.get_dummies(m.lineage_id,dtype=float).to_numpy();S=pd.get_dummies(m.dataset_id,dtype=float).to_numpy()[:,1:];Q=m[['completeness','contamination']].to_numpy();Q=(Q-Q.mean(0))/Q.std(0,ddof=1);xr=np.column_stack([Z,Q]);x=np.column_stack([xr,S]);rank=xr.shape[1];assert np.linalg.matrix_rank(x)==x.shape[1]
groups=m.groupby('lineage_id',sort=True).indices;indices=list(groups.values());rng=np.random.default_rng(20261002);w=rng.exponential(size=(2000,len(indices)));np.save(O/'strict_BB_weights.npy',w)
for branch in ['core65','restored74']:
 y=b[b.branch==branch].pivot(index='tip_id',columns='module',values='breadth').loc[m.index,cols].to_numpy();xx=np.array([x[i].T@x[i] for i in indices]);xy=np.array([x[i].T@y[i] for i in indices]);yy=np.array([(y[i]**2).sum(0) for i in indices]);ys=np.array([y[i].sum(0) for i in indices]);ns=np.array([len(i) for i in indices])
 def fit(weights):
  XX=np.einsum('bg,gij->bij',weights,xx);XY=np.einsum('bg,gik->bik',weights,xy);YY=weights@yy;YS=weights@ys;N=weights@ns;var=YY/N[:,None]-(YS/N[:,None])**2;assert np.all(var>0)
  bf=np.linalg.solve(XX,XY);br=np.linalg.solve(XX[:,:rank,:rank],XY[:,:rank]);sf=(YY-np.einsum('bik,bik->bk',bf,XY))/var;sr=(YY-np.einsum('bik,bik->bk',br,XY[:,:rank]))/var;assert np.all(sf>=-1e-6) and np.all(sr>=sf-1e-6)
  return sr,sf
 sr,sf=fit(np.ones((1,len(indices))));ma=1-sf[0,:5].sum()/sr[0,:5].sum();en=1-sf[0,5:].sum()/sr[0,5:].sum();z=(y-y.mean(0))/y.std(0,ddof=1);r0=z-xr@np.linalg.lstsq(xr,z,rcond=None)[0];r1=z-x@np.linalg.lstsq(x,z,rcond=None)[0];qr=np.array([1-(r1[:,:5]**2).sum()/(r0[:,:5]**2).sum(),1-(r1[:,5:]**2).sum()/(r0[:,5:]**2).sum()]);err=float(np.max(abs(qr-[ma,en])));assert err<1e-9
 draws=[]
 for start in range(0,2000,100):
  a,c=fit(w[start:start+100]);draws.extend((1-c[:,:5].sum(1)/a[:,:5].sum(1))-(1-c[:,5:].sum(1)/a[:,5:].sum(1)))
 d=np.array(draws);np.save(O/f'{branch}_strict_delta_BB.npy',d);results.append(dict(branch=branch,block_kind=args.block_kind,n=len(m),lineage_blocks=len(indices),catalogues=m.dataset_id.nunique(),maintenance_partial_R2=ma,energy_partial_R2=en,delta_point=ma-en,delta_BB_mean=float(d.mean()),delta_q025=float(np.quantile(d,.025)),delta_q975=float(np.quantile(d,.975)),Pr_delta_positive=float((d>0).mean()),draws=2000,scope='Paired block Bayesian bootstrap conditional on observed catalogues; not environmental-effect posterior'))
 checks.append(dict(branch=branch,reduced_rank=rank,full_rank=x.shape[1],QR_point_max_error=err))
 for j,col in enumerate(cols):decomp.append(dict(branch=branch,module=col,source_partial_R2=float(1-sf[0,j]/sr[0,j]),mean_breadth=float(y[:,j].mean()),sd_breadth=float(y[:,j].std(ddof=1))))
pd.DataFrame(results).to_csv(O/'strict_paired_uncertainty.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(O/'strict_model_validation.tsv',sep='\t',index=False);pd.DataFrame(decomp).to_csv(O/'strict_module_decomposition.tsv',sep='\t',index=False);print(pd.DataFrame(results).to_string(index=False))
