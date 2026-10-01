"""Exploratory localization-aware refinement; no HMM score is a function probability.
Greedy best i-Evalue domains, retaining near-score alternatives separately.
Not an exact reproduction of run_dbcan's adjacent-only overlap filter.
"""
from pathlib import Path
import pandas as pd,json
import argparse
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
H=args.output;H.mkdir(parents=True,exist_ok=True)
INPUT=Path(__file__).resolve().parent/'inputs/subfamily_passing_domains.tsv.gz'
x=pd.read_csv(INPUT,sep='\t',keep_default_na=False,low_memory=False)
rows=[];selected=[]
def overlap(a,b):
 n=max(0,min(a['ali_to'],b['ali_to'])-max(a['ali_from'],b['ali_from'])+1)
 return n/min(a['ali_to']-a['ali_from']+1,b['ali_to']-b['ali_from']+1)>.5
for protein,g in x.groupby('protein',sort=True):
 records=g.sort_values(['domain_e','domain_score','subfamily','ali_from'],ascending=[True,False,True,True]).to_dict('records');keepers=[]
 for r in records:
  if not any(overlap(r,k) for k in keepers):keepers.append(r)
 for i,k in enumerate(keepers,1):
  selected.append(dict(k,domain_index=i))
  candidates=[r for r in records if overlap(r,k)]
  for margin in [0,10,20]:
   # E-value selects anchor; score tolerance is sensitivity, not calibrated confidence.
   options=[r for r in candidates if r['domain_score']>=k['domain_score']-margin]
   if k not in options:options.append(k)
   labels=sorted({s for r in options for s in r['substrates'].split(';') if s});unknown=any(not r['substrates'] for r in options)
   rows.append(dict(genome=k['genome'],protein=protein,domain_index=i,anchor_family=k['family'],anchor_subfamily=k['subfamily'],anchor_score=k['domain_score'],score_margin=margin,alternative_models=len({r['subfamily'] for r in options}),has_unknown=unknown,labels=';'.join(labels),status='unmapped' if not labels else ('single_label_candidate' if len(labels)==1 and not unknown else 'ambiguous'),HMO_candidate='human milk polysaccharide' in labels))
a=pd.DataFrame(selected);a.to_csv(H/'nonoverlap_anchor_domains.tsv',sep='\t',index=False)
r=pd.DataFrame(rows);r.to_csv(H/'domain_competition_sensitivity.tsv',sep='\t',index=False)
q=[]
for margin,g in r.groupby('score_margin'):
 q.append(dict(score_margin=int(margin),domains=len(g),proteins=g.protein.nunique(),status=g.status.value_counts().to_dict(),single_label_non_CBM_GT=int(((g.status=='single_label_candidate')&~g.anchor_family.str.startswith(('CBM','GT'))).sum()),HMO_domains=int(g.HMO_candidate.sum()),HMO_genomes=g[g.HMO_candidate].genome.nunique()))
(H/'domain_competition_summary.json').write_text(json.dumps(q,indent=2));print(json.dumps(q,indent=2))
