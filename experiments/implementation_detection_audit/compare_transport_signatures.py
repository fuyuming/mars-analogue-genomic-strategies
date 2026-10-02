from pathlib import Path
import pandas as pd,json,numpy as np
H=Path(__file__).resolve().parent;O=H/'host';x=pd.read_csv(O/'all5477_expanded10_presence.tsv.gz',sep='\t',index_col=0);m=pd.read_csv(O/'signature_calls.tsv',sep='\t');rows=[];fc=[];allrows=[]
for q,d in m.groupby('quality'):
 d=d.copy();d['OpuABC']=x.loc[d.genome_id,['K05845','K05846','K05847']].eq(1).all(axis=1).to_numpy();d['BCCT']=x.loc[d.genome_id,'K02168'].eq(1).to_numpy();d['OpuD_BetL']=x.loc[d.genome_id,'K05020'].eq(1).to_numpy();d['pair_state']=np.select([d.ProVWX&d.OpuABC,d.ProVWX&~d.OpuABC,~d.ProVWX&d.OpuABC],['both','ProU_signature_only','Opu_signature_only'],'neither');allrows.append(d)
 for domain,g in d.groupby('expected_domain'):
  for s,c in g.pair_state.value_counts().items():rows.append(dict(quality=q,domain=domain,state=s,n=c))
  for fam,g2 in g[g.family.ne('Unclassified')].groupby('family'):
   a=g2[g2.pair_state.eq('ProU_signature_only')];b=g2[g2.pair_state.eq('Opu_signature_only')];ok=len(a)>=3 and len(b)>=3 and a.dataset_id.nunique()>=2 and b.dataset_id.nunique()>=2
   fc.append(dict(quality=q,domain=domain,family=fam,n=len(g2),ProU_only=len(a),Opu_only=len(b),ProU_catalogues=a.dataset_id.nunique(),Opu_catalogues=b.dataset_id.nunique(),both=int(g2.pair_state.eq('both').sum()),coarse_comparability=ok))
a=pd.concat(allrows);a.to_csv(O/'all_final_transport_calls.tsv',sep='\t',index=False);pd.DataFrame(rows).to_csv(O/'transport_pair_states.tsv',sep='\t',index=False);f=pd.DataFrame(fc);f.to_csv(O/'transport_family_comparability.tsv',sep='\t',index=False);f[f.coarse_comparability].to_csv(O/'transport_comparable_families.tsv',sep='\t',index=False)
h=pd.DataFrame(json.loads((O/'additional_markers_new556.json').read_text())['protein_hits']);shared=[];P={'K02000','K02001','K02002'};U={'K05845','K05846','K05847'}
for (g,p),v in h.groupby(['genome_id','protein']):
 ks=set(v.KO)
 if ks&P and ks&U:shared.append(dict(genome_id=g,protein=p,ProU_KOs=';'.join(sorted(ks&P)),Opu_KOs=';'.join(sorted(ks&U))))
pd.DataFrame(shared,columns=['genome_id','protein','ProU_KOs','Opu_KOs']).to_csv(O/'shared_protein_route_assignments.tsv',sep='\t',index=False)
print(pd.DataFrame(rows).to_string(index=False));print(f[f.coarse_comparability].to_string(index=False));print('shared protein assignments',len(shared),'genomes',len({r['genome_id'] for r in shared}))
