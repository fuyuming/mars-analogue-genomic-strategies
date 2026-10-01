from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host';N=R/'48_MAG_portfolio_inference_20261001/host/new_catalogue'
RULES={'KdpABC':['K01546','K01547','K01548'],'KdpDE':['K07646','K07667'],'ProVWX':['K02000','K02001','K02002'],'EctABC':['K06718','K00836','K06720'],'EctD':['K10674'],'BetAB':['K00108','K00130']}
old=pd.read_csv(A/'restored74_presence.csv.gz',usecols=['genome_id','KO','present']);old=old[old.KO.isin(sum(RULES.values(),[]))].pivot(index='genome_id',columns='KO',values='present').astype(int)
new=pd.read_csv(N/'restored74_presence.tsv.gz',sep='\t',index_col=0);meta=pd.read_csv(N/'desert826_endpoint_metadata.tsv',sep='\t').set_index('genome_id')
site_union=meta.groupby('ani_cluster_95').site_labels.agg(lambda a:';'.join(sorted(set(a))));site_n=meta.groupby('ani_cluster_95').site_labels.nunique();meta['ANI_site_union']=meta.ani_cluster_95.map(site_union);meta['ANI_n_sites']=meta.ani_cluster_95.map(site_n)
meta['selection_score']=meta.checkm2_completeness-5*meta.checkm2_contamination
cohorts={};audit=[]
for strict in [False,True]:
 st='strict' if strict else 'primary'
 om=pd.read_csv(R/f'48_MAG_portfolio_inference_20261001/host/restored74_{st}_cohort.tsv',sep='\t').set_index('genome_id');cohorts['original_'+st]=(old.loc[om.index],om,'lineage_id')
 eligible=(meta.domain=='Bacteria') & (~meta.background_flagged) & (~meta.shares_ANI95_with_any_existing)
 if strict:eligible &= meta.passes_strict_90_5
 em=meta[eligible].copy();em['id_sort']=em.index;em=em.sort_values(['selection_score','contig_n50','id_sort'],ascending=[False,False,True]).drop_duplicates('ani_cluster_95');assert em.ani_cluster_95.nunique()==len(em)
 em.to_csv(O/f'new_desert_{st}_independent_representatives.tsv',sep='\t')
 cohorts['new_desert_'+st]=(new.loc[em.index],em,'site_labels')
 audit.append({'cohort':'new_desert_'+st,'n':len(em),'sites':em.site_labels.nunique(),'classes':em['class'].nunique(),'multiple_site_ANI_clusters':int((em.ANI_n_sites>1).sum()),'selection':'Bacteria; no flagged genome; no ANI95 overlap with any existing4095; one rep per ANI95 quality score then N50 then ID; quality filter before representative choice'})
full=[];summ=[];pairs=[]
for cohort,(y,m,block) in cohorts.items():
 statuses=pd.DataFrame(index=y.index)
 for system,ks in RULES.items():
  hits=y[ks].sum(axis=1);state=np.select([hits==len(ks),hits==0],['panel_complete','not_detected'],default='partial');statuses[system]=state
  nr={'cohort':cohort,'system':system,'n':len(y),'not_detected':int((hits==0).sum()),'partial':int(((hits>0)&(hits<len(ks))).sum()),'panel_complete':int((hits==len(ks)).sum()),'n_KO':len(ks),'complete_given_any':float((hits==len(ks)).sum()/max(1,(hits>0).sum()))}
  groups=m.groupby(block).groups;counts=[]
  for group,ix in groups.items():
   z=hits.loc[ix];counts.append([int((z==len(ks)).sum()),int((z>0).sum())]);full.append({'cohort':cohort,'system':system,'group':group,'n':len(z),'panel_complete':int((z==len(ks)).sum()),'any_hit':int((z>0).sum())})
  counts=np.array(counts);rng=np.random.default_rng(20261002);w=rng.exponential(size=(2000,len(counts)));draw=(w@counts[:,0])/(w@counts[:,1]);nr['conditional_reweight_q025'],nr['conditional_reweight_q975']=np.quantile(draw,[.025,.975]);nr['reweight_unit']=block;summ.append(nr)
 statuses.join(m).to_csv(O/f'{cohort}_system_states.tsv',sep='\t')
 for a,b in [('KdpABC','KdpDE'),('EctABC','EctD'),('EctABC','ProVWX'),('BetAB','ProVWX')]:
  for sa in ['not_detected','partial','panel_complete']:
   for sb in ['not_detected','partial','panel_complete']:
    ix=(statuses[a]==sa)&(statuses[b]==sb);pairs.append({'cohort':cohort,'systemA':a,'systemB':b,'stateA':sa,'stateB':sb,'n_MAG':int(ix.sum()),'n_groups':int(m.loc[ix,block].nunique())})
pd.DataFrame(summ).to_csv(O/'system_state_summary.tsv',sep='\t',index=False);pd.DataFrame(full).to_csv(O/'system_state_group_counts.tsv',sep='\t',index=False);pd.DataFrame(pairs).to_csv(O/'paired_system_configurations.tsv',sep='\t',index=False);pd.DataFrame(audit).to_csv(O/'external_selection_audit.tsv',sep='\t',index=False)
(O/'system_rules_provisional.json').write_text(json.dumps({'rules':RULES,'meaning':'KO co-detection only; no adjacency, completeness or activity claim; functional interpretation pending independent definitions'},indent=2))
print(pd.DataFrame(audit).to_string(index=False));print(pd.DataFrame(summ)[['cohort','system','n','not_detected','partial','panel_complete','complete_given_any']].to_string(index=False))
