from pathlib import Path
from itertools import product
import pandas as pd,numpy as np,json,re,hashlib
H=Path(__file__).resolve().parent;R=H.parent;O=H/'host';P=R/'49_system_configuration_20261001/host'
rules=json.loads((P/'system_rules_provisional.json').read_text())['rules']
old=pd.read_csv(O/'original_target_proteins.tsv',sep='\t')
parts=old.protein_id.str.extract(r'^(.*)_(\d+)$');old['contig']=parts[0];old['cds_rank']=pd.to_numeric(parts[1]);assert old.contig.notna().all()
new=pd.read_csv(O/'desert_target_coordinates.tsv',sep='\t',dtype={'partial':str});new['genome_id']='desert2025__'+new.genome_id
allrows=[];examples=[];summ=[];bygroup=[];patterns=[]
for cohort in ['original_primary','original_strict','new_desert_primary','new_desert_strict']:
 m=pd.read_csv(P/f'{cohort}_system_states.tsv',sep='\t').set_index('genome_id');isnew=cohort.startswith('new');data=(new if isnew else old);groups={g:d for g,d in data[data.genome_id.isin(m.index)].groupby('genome_id')}
 for gid,mr in m.iterrows():
  d=groups.get(gid,data.iloc[:0]);unit=mr.site_labels if isnew else mr.lineage_id;source=mr.site_labels if isnew else mr.dataset_id
  for system,ks in rules.items():
   hit=d[d.KO.isin(ks)].copy();obs=set(hit.KO);codetected=set(ks)<=obs
   assert codetected==(mr[system]=='panel_complete'),(cohort,gid,system,'co-detection mismatch')
   assert bool(obs)==(mr[system]!='not_detected'),(cohort,gid,system,'any-hit mismatch')
   row=dict(cohort=cohort,genome_id=gid,source_or_site=source,unit=unit,system=system,n_components=len(ks),any_hit=bool(obs),codetected=codetected,same_contig=False,order_compact=False,compact_5kb=False if isnew else None,compact_10kb=False if isnew else None,compact_20kb=False if isnew else None,min_span_bp=np.nan,min_span_cds=np.nan,any_target_partial=bool((hit.partial!='00').any()) if isnew else None)
   if codetected and len(ks)>1:
    candidates=[]
    for contig,g in hit.groupby('contig'):
     if not set(ks)<=set(g.KO):continue
     for combo in product(*(g[g.KO==k].to_dict('records') for k in ks)):
      if len({x['protein_id'] for x in combo})<len(ks):continue
      rankspan=max(x['cds_rank'] for x in combo)-min(x['cds_rank'] for x in combo)+1
      span=max(x['end'] for x in combo)-min(x['start'] for x in combo)+1 if isnew else np.nan
      ss=len({x['strand'] for x in combo})==1 if isnew else None
      candidates.append((rankspan,span,ss,contig,combo))
    if candidates:
     row['same_contig']=True;row['min_span_cds']=min(x[0] for x in candidates);row['order_compact']=row['min_span_cds']<=len(ks)+2
     if isnew:
      row['min_span_bp']=min(x[1] for x in candidates)
      for th in [5,10,20]:row[f'compact_{th}kb']=any(x[2] and x[1]<=th*1000 for x in candidates)
      best=sorted(candidates,key=lambda x:(not x[2],x[1],x[3]))[0]
      if cohort=='new_desert_strict' and row['compact_10kb']:
       for x in best[4]:examples.append({**x,'system':system,'source_or_site':source,'span_bp':best[1]})
   allrows.append(row)
   pat=''.join('1' if k in obs else '0' for k in ks)
   patterns.append(dict(cohort=cohort,genome_id=gid,system=system,pattern=pat,KO_order=';'.join(ks),source_or_site=source))
out=pd.DataFrame(allrows);out.to_csv(O/'system_context_by_MAG.tsv',sep='\t',index=False)
for (cohort,system),g in out.groupby(['cohort','system'],sort=False):
 if system=='EctD':continue
 v=g[g.codetected];counts=v.groupby('unit')[['same_contig','order_compact','compact_5kb','compact_10kb','compact_20kb']].sum();den=v.groupby('unit').size()
 for metric in counts:
  if cohort.startswith('original') and metric.startswith('compact_'):continue
  n=int(v[metric].sum());total=len(v);lo=hi=np.nan
  if total:
   rng=np.random.default_rng(20261003);w=rng.exponential(size=(2000,len(counts)));draw=(w@counts[metric].values)/(w@den.values);lo,hi=np.quantile(draw,[.025,.975])
  summ.append(dict(cohort=cohort,system=system,metric=metric,n=n,denominator=total,fraction=n/total if total else np.nan,conditional_q025=lo,conditional_q975=hi,group_count=len(counts),interpretation='verified coordinates' if cohort.startswith('new') else 'protein-ID contig/order proxy only'))
  for unit in counts.index:bygroup.append(dict(cohort=cohort,system=system,metric=metric,unit=unit,n=int(counts.loc[unit,metric]),denominator=int(den.loc[unit])))
pd.DataFrame(summ).to_csv(O/'system_context_summary.tsv',sep='\t',index=False)
pd.DataFrame(bygroup).to_csv(O/'system_context_group_counts.tsv',sep='\t',index=False)
pd.DataFrame(examples).to_csv(O/'strict_compact_locus_candidates.tsv',sep='\t',index=False)
pd.DataFrame(patterns).groupby(['cohort','system','pattern','KO_order']).size().rename('n').reset_index().to_csv(O/'system_patterns.tsv',sep='\t',index=False)
audit={'all_four_cohort_marker_states_match_round49':True,'original_ID_parseable':int(old.contig.notna().sum()),'original_target_rows':len(old),'desert_target_rows':len(new),'all_desert_targets_mapped_to_GFF':True,'no_reannotation':True,'original_coordinates_not_verified':True,'scope':'post-result molecular-context diagnostic; no new ecology model','inputs':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [O/'original_target_proteins.tsv',O/'desert_target_coordinates.tsv',P/'system_rules_provisional.json']}}
(O/'context_validation.json').write_text(json.dumps(audit,indent=2))
print(pd.DataFrame(summ).query("metric in ['same_contig','order_compact','compact_10kb']")[['cohort','system','metric','n','denominator']].to_string(index=False))
