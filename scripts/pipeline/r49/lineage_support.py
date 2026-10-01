from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pandas as pd,numpy as np,json
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host';N=R/'48_MAG_portfolio_inference_20261001/host/new_catalogue'
RULES=json.loads((O/'system_rules_provisional.json').read_text())['rules'];panel=pd.read_csv(A/'restored74_panel.tsv',sep='\t');panel=panel[panel.module=='osmotic_desiccation_salt'];raw=pd.read_csv(A/'restored74_presence.csv.gz');tax=raw[['genome_id','class','order','family']].drop_duplicates().set_index('genome_id');orig=raw[raw.KO.isin(panel.KO)].pivot(index='genome_id',columns='KO',values='present').astype(float)
new=pd.read_csv(N/'restored74_presence.tsv.gz',sep='\t',index_col=0).astype(float)
allstats=[];support=[];groupmeans=[];patterns=[]
for cohort in ['original_primary','original_strict','new_desert_primary','new_desert_strict']:
 m=pd.read_csv(O/f'{cohort}_system_states.tsv',sep='\t').set_index('genome_id');isold=cohort.startswith('original');source='dataset_id' if isold else 'site_labels';y=(orig if isold else new).loc[m.index]
 if isold:m=m.join(tax)
 response={}
 for sub,pp in panel.groupby('submodule'):response[sub]=y[list(pp.KO)].mean(axis=1)
 for name,ks in RULES.items():
  response[name+'_breadth']=y[ks].mean(axis=1);response[name+'_codetected']=y[ks].min(axis=1)
 ys=pd.DataFrame(response)
 for src,ix in m.groupby(source).groups.items():
  for col in ys:groupmeans.append({'cohort':cohort,'source_or_site':src,'endpoint':col,'n':len(ix),'mean':ys.loc[ix,col].mean()})
 for name,ks in RULES.items():
  v=y[ks].astype(int).astype(str).agg(''.join,axis=1)
  for (src,pat),g in pd.DataFrame({'source':m[source],'pattern':v}).groupby(['source','pattern']):patterns.append({'cohort':cohort,'system':name,'KO_order':';'.join(ks),'source_or_site':src,'pattern':pat,'n':len(g)})
 ranks=['lineage_id','class','order','family'] if isold else ['class','order','family']
 for rank in ranks:
  counts=m.groupby([rank,source]).size().unstack(fill_value=0);keep=counts.index[(counts>=5).sum(axis=1)>=2]
  # All representatives in an eligible lineage included; threshold establishes common support only.
  take=m[rank].isin(keep)&m[rank].notna()&~m[rank].astype(str).str.lower().str.startswith('unclassified');mm=m[take];z=ys.loc[mm.index].to_numpy(float)
  support.append({'cohort':cohort,'rank':rank,'n':len(mm),'lineages':mm[rank].nunique(),'source_or_sites':mm[source].nunique(),'gate':'at least2 groups with>=5 genomes in each; all members of eligible lineage retained'})
  if len(mm)<20 or mm[source].nunique()<2:continue
  q=mm[['checkm2_completeness','checkm2_contamination']].to_numpy();q=(q-q.mean(0))/q.std(0);xr=np.column_stack([pd.get_dummies(mm[rank],dtype=float),q]);xf=np.column_stack([xr,pd.get_dummies(mm[source],dtype=float).iloc[:,1:]]);rr=np.linalg.matrix_rank(xr);rf=np.linalg.matrix_rank(xf);er=z-xr@np.linalg.lstsq(xr,z,rcond=None)[0];ef=z-xf@np.linalg.lstsq(xf,z,rcond=None)[0];sr=(er**2).sum(0);sf=(ef**2).sum(0)
  for j,col in enumerate(ys):
   valid=sr[j]>1e-9
   allstats.append({'cohort':cohort,'rank':rank,'endpoint':col,'n':len(mm),'lineages':mm[rank].nunique(),'groups':mm[source].nunique(),'rank_reduced':rr,'rank_full':rf,'positive':int((z[:,j]>0).sum()),'positive_groups':int(mm.loc[z[:,j]>0,source].nunique()),'zero_groups':int(mm.loc[z[:,j]==0,source].nunique()),'limited_observed_state_support':bool((z[:,j]>0).sum()<10 or mm.loc[z[:,j]>0,source].nunique()<3 or (z[:,j]==0).sum()<10),'source_or_site_partial_r2':1-sf[j]/sr[j] if valid else np.nan,'df_adjusted_partial_r2':1-(sf[j]/(len(mm)-rf))/(sr[j]/(len(mm)-rr)) if valid and len(mm)>rf else np.nan,'status':'descriptive_group_projection' if valid else 'invariant_after_adjustment'})
pd.DataFrame(allstats).to_csv(O/'lineage_group_associations.tsv',sep='\t',index=False);pd.DataFrame(support).to_csv(O/'lineage_group_support.tsv',sep='\t',index=False);pd.DataFrame(groupmeans).to_csv(O/'source_site_system_means.tsv',sep='\t',index=False);pd.DataFrame(patterns).to_csv(O/'system_marker_patterns_by_source.tsv',sep='\t',index=False)
print(pd.DataFrame(support).to_string(index=False));print(pd.DataFrame(allstats).query("endpoint in ['KdpABC_breadth','KdpDE_breadth','potassium_homeostasis']")[['cohort','rank','endpoint','n','source_or_site_partial_r2','df_adjusted_partial_r2']].to_string(index=False))
