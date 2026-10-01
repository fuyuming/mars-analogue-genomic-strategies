import analyse_portfolios as a
import numpy as np,pandas as pd
O=a.OUT
for branch in ['core65','restored74']:
 pan=pd.read_csv(a.A/f'{branch}_panel.tsv',sep='\t');pan=pan[pan.module=='osmotic_desiccation_salt']
 ko=pd.read_csv(a.A/f'{branch}_presence.csv.gz',usecols=['genome_id','KO','present']);ko=ko[ko.KO.isin(pan.KO)]
 for strict in [False,True]:
  meta,y,x,r,C,stats=a.prepare(branch,strict);y=y.copy()
  ky=ko.pivot(index='genome_id',columns='KO',values='present').loc[meta.index].astype(float)
  y[:,3]=np.column_stack([ky[list(g.KO)].mean(axis=1) for _,g in pan.groupby('submodule')]).mean(1)
  cells=meta.reset_index().groupby(['lineage_id','dataset_id'],sort=True).indices
  stats[1]=np.array([x[cells[c]].T@y[cells[c]] for c in C]);stats[2]=np.array([(y[cells[c]]**2).sum(0) for c in C]);stats[3]=np.array([y[cells[c]].sum(0) for c in C])
  li=sorted(meta.lineage_id.unique());si=sorted(meta.dataset_id.unique());lc=np.array([li.index(c[0]) for c in C]);sc=np.array([si.index(c[1]) for c in C]);rng=np.random.default_rng(a.SEED);lw=rng.exponential(size=(a.B,len(li)));sw=rng.exponential(size=(a.B,len(si)))
  for scheme,w in [('lineage_BB',lw[:,lc]),('lineage_source_product_sensitivity',lw[:,lc]*sw[:,sc])]:
   chunks=[a.fit_weights(w[i:i+100],stats,r) for i in range(0,a.B,100)];sr=np.concatenate([c[0] for c in chunks]);sf=np.concatenate([c[1] for c in chunks]);a.summary(branch,strict,scheme,sr,sf,'equal_osmotic_submodules',list(range(5)),[5,6]);np.savez_compressed(O/f'{branch}_{strict}_{scheme}_equal_osmotic_draws.npz',sr=sr,sf=sf)
pd.DataFrame(a.results).to_csv(O/'equal_submodule_uncertainty.tsv',sep='\t',index=False)
print(pd.DataFrame(a.results).to_string(index=False))
