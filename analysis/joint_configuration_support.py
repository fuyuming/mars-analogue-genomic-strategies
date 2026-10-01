"""Descriptive support for genomic marker combinations; no activity/causality inference."""
import argparse,json
from pathlib import Path
import pandas as pd
RULES={'KdpABC':['K01546','K01547','K01548'],'ProVWX':['K02000','K02001','K02002'],'EctABC':['K06718','K00836','K06720'],'BetAB':['K00108','K00130'],'CoxCut':['K03518','K03519','K03520']}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--revision-root',type=Path,required=True);ap.add_argument('--original-presence',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True);a=ap.parse_args();a.outdir.mkdir(parents=True,exist_ok=True)
 r=a.revision_root;h=r/'49_system_configuration_20261001/host';raw=pd.read_csv(a.original_presence);tax=raw[['genome_id','family','genus']].drop_duplicates().set_index('genome_id');assert not tax.index.duplicated().any()
 old=raw.pivot(index='genome_id',columns='KO',values='present');new=pd.read_csv(r/'48_MAG_portfolio_inference_20261001/host/new_catalogue/restored74_presence.tsv.gz',sep='\t',index_col=0)
 summaries=[];cells=[];support=[];allstates=[]
 for cohort in ['original_primary','original_strict','new_desert_primary','new_desert_strict']:
  states=pd.read_csv(h/f'{cohort}_system_states.tsv',sep='\t').set_index('genome_id')
  if cohort.startswith('original'):
   meta=states.join(tax);meta['unit']=meta.dataset_id;mat=old.loc[meta.index]
  else:
   suffix=cohort.replace('new_desert_','');meta=pd.read_csv(h/f'new_desert_{suffix}_independent_representatives.tsv',sep='\t').set_index('genome_id');meta=meta.loc[states.index];meta['unit']=meta.site_labels;mat=new.loc[meta.index]
  assert mat.index.is_unique and mat.loc[:,sum(RULES.values(),[])].notna().all().all()
  fam=meta.family.fillna('');known=~fam.isin(['','f__','unclassified','unassigned','nan']);fam=fam.where(known)
  calls={}
  for name,kos in RULES.items():
   hits=mat[kos].sum(axis=1);calls[name]=hits.eq(len(kos));c=hits.map(lambda n:'all_detected' if n==len(kos) else ('partial' if n>0 else 'not_detected'))
   if name!='CoxCut':assert (calls[name]==states[name].eq('panel_complete')).all()
   for gid in meta.index:allstates.append({'cohort':cohort,'genome_id':gid,'family':fam.loc[gid],'unit':meta.loc[gid,'unit'],'system':name,'state':c.loc[gid]})
  for name in list(RULES)[:-1]:
   d=pd.DataFrame({'maintenance':calls[name].astype(int),'energy':calls['CoxCut'].astype(int),'family':fam,'unit':meta.unit});d['joint']=2*d.maintenance+d.energy
   ct=d.joint.value_counts();eligible=[];familyexp=0;familyobs=0;nk=0
   for family,g in d.dropna(subset=['family']).groupby('family'):
    counts=[int((g.joint==j).sum()) for j in range(4)];units=g.unit.value_counts();n=len(g);exp=float(g.maintenance.sum()*g.energy.sum()/n)
    familyexp+=exp;familyobs+=counts[3];nk+=n
    cells.append({'cohort':cohort,'pair':name+'__CoxCut','family':family,'n':n,'n00':counts[0],'n01':counts[1],'n10':counts[2],'n11':counts[3],'source_or_site_groups':len(units),'groups_ge3':int((units>=3).sum()),'expected_both_within_family':exp})
    for gate in [2,3,5]:
     if min(counts)>=gate and (units>=3).sum()>=2:eligible.append((gate,family,n,len(units)))
   for gate in [2,3,5]:
    ee=[x for x in eligible if x[0]==gate];support.append({'cohort':cohort,'pair':name+'__CoxCut','min_per_joint_state':gate,'eligible_families':len(ee),'eligible_genomes':sum(x[2] for x in ee),'families':';'.join(x[1] for x in ee)})
   summaries.append({'cohort':cohort,'pair':name+'__CoxCut','n':len(d),'n00':int(ct.get(0,0)),'n01':int(ct.get(1,0)),'n10':int(ct.get(2,0)),'n11':int(ct.get(3,0)),'expected_both_pooled':float(d.maintenance.sum()*d.energy.sum()/len(d)),'family_assigned_n':nk,'both_family_assigned':familyobs,'expected_both_family_conditioned':familyexp,'expected_both_pooled_family_assigned':float(d.loc[known,'maintenance'].sum()*d.loc[known,'energy'].sum()/max(1,nk)),'unassigned_family_n':len(d)-nk})
 for name,rows in [('joint_summary',summaries),('family_joint_cells',cells),('joint_support',support),('genome_system_states',allstates)]:pd.DataFrame(rows).to_csv(a.outdir/(name+'.tsv'),sep='\t',index=False)
 (a.outdir/'validation.json').write_text(json.dumps({'all_prior_maintenance_calls_match':True,'pairs':4,'cohorts':4,'CO_activity_claimed':False,'environment_model_fitted':False,'thresholds':[2,3,5]},indent=2));print(pd.DataFrame(support).query('min_per_joint_state==3').to_string(index=False))
if __name__=='__main__':main()
