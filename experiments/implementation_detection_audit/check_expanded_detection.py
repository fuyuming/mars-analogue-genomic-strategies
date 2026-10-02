from pathlib import Path
import pandas as pd,numpy as np,json
H=Path(__file__).resolve().parent;O=H/'host';r=json.loads((O/'additional_markers_new556.json').read_text());x=pd.DataFrame(r['genomes']).set_index('genome_id');calls=pd.read_csv(O/'signature_calls.tsv',sep='\t');res=pd.read_csv(O/'source_checked_resource_signature_calls.tsv',sep='\t');out=[];allrows=[];rr=[]
sigs={'BCCT_K02168':['K02168'],'OpuD_BetL_K05020':['K05020'],'Opu_ABC_K05845_6_7':['K05845','K05846','K05847']}
for q,m in calls.groupby('quality'):
 d=m[m.genome_id.isin(x.index)].copy().set_index('genome_id');z=x.loc[d.index];assert len(d)==(499 if q=='primary' else 123)
 for key,ks in sigs.items():d[key]=z[ks].eq(1).all(axis=1)
 d['any_additional_signature']=d[list(sigs)].any(axis=1);d['old_or_additional']=d.ProVWX|d.any_additional_signature
 for (source,domain),g in d.groupby(['dataset_id','expected_domain']):
  out.append(dict(quality=q,dataset_id=source,domain=domain,n=len(g),ProVWX_all=int(g.ProVWX.sum()),**{k:int(g[k].sum()) for k in sigs},any_additional=int(g.any_additional_signature.sum()),added_when_ProVWX_not_all=int((g.any_additional_signature&~g.ProVWX).sum()),union_signatures=int(g.old_or_additional.sum())))
 joined=d.reset_index().merge(res[res.quality.eq(q)][['genome_id','sample_id','region_group']],on='genome_id',how='inner',validate='one_to_one')
 for key in ['ProVWX',*sigs,'old_or_additional']:
  rr.append(dict(quality=q,signature=key,n=len(joined),detected=int(joined[key].sum()),samples_with=int(joined.loc[joined[key],'sample_id'].nunique()),samples_without=int(joined.loc[~joined[key],'sample_id'].nunique()),families_with_both_states=int(joined.groupby('family')[key].nunique().eq(2).sum())))
 allrows.append(d.reset_index());joined.to_csv(O/f'{q}_expanded_resource_calls.tsv',sep='\t',index=False)
pd.concat(allrows).to_csv(O/'new_cohort_expanded_signatures.tsv',sep='\t',index=False);pd.DataFrame(out).to_csv(O/'expanded_detection_summary.tsv',sep='\t',index=False);pd.DataFrame(rr).to_csv(O/'expanded_resource_variation.tsv',sep='\t',index=False);print(pd.DataFrame(out).to_string(index=False));print(pd.DataFrame(rr).to_string(index=False))
