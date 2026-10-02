from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent;O=H/'host';d=pd.read_csv(O/'all_final_transport_calls.tsv',sep='\t');rows=[]
for (q,domain,genus),g in d[d.genus.ne('Unclassified')].groupby(['quality','expected_domain','genus']):
 a=g[g.pair_state.eq('ProU_signature_only')];b=g[g.pair_state.eq('Opu_signature_only')]
 if len(a) and len(b):rows.append(dict(quality=q,domain=domain,genus=genus,n=len(g),ProU_only=len(a),Opu_only=len(b),ProU_catalogues=a.dataset_id.nunique(),Opu_catalogues=b.dataset_id.nunique(),coarse_comparability=len(a)>=3 and len(b)>=3 and a.dataset_id.nunique()>=2 and b.dataset_id.nunique()>=2))
r=pd.DataFrame(rows);r.to_csv(O/'transport_genus_both_states.tsv',sep='\t',index=False);print(r[r.coarse_comparability].to_string(index=False));print('Strict both-states genera:',r[r.quality.eq('strict')].to_string(index=False))
