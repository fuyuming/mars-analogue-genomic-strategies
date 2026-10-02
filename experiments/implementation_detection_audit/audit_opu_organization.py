from pathlib import Path
import json,pandas as pd,itertools,re
H=Path(__file__).resolve().parent;O=H/'host';raw=json.loads((O/'additional_markers_new556.json').read_text());coords=pd.DataFrame(json.loads((O/'opu_original_coordinates.json').read_text()));hits=pd.DataFrame(raw['protein_hits']);d=hits[hits.KO.isin(['K05845','K05846','K05847'])].merge(coords,on=['genome_id','protein'],how='inner',validate='many_to_one');d['partial']=d.header.str.extract(r'partial=(\d+)');d.to_csv(O/'opu_protein_evidence.tsv',sep='\t',index=False);out=[];triples=[]
for genome,g in d.groupby('genome_id'):
 cand=[]
 for (contig,strand),a in g.groupby(['contig','strand']):
  if set(a.KO)!={'K05845','K05846','K05847'}:continue
  for combo in itertools.product(*[list(a[a.KO.eq(k)].to_dict('records')) for k in ['K05845','K05846','K05847']]):
   if len({z['protein'] for z in combo})<3:continue
   ordered=sorted(combo,key=lambda z:z['start']);span=max(z['orf_index'] for z in combo)-min(z['orf_index'] for z in combo);gap=max(max(0,ordered[i+1]['start']-ordered[i]['end']-1) for i in [0,1]);v=dict(genome_id=genome,contig=contig,strand=int(strand),orf_span=int(span),max_gap_bp=int(gap),all_partial00=all(z['partial']=='00' for z in combo),within_declared_window=span<=5 and gap<=2000,**{z['KO']+'_protein':z['protein'] for z in combo});cand.append(v)
 cand.sort(key=lambda c:(not c['within_declared_window'],not c['all_partial00'],c['orf_span'],c['max_gap_bp'],c['contig']))
 rec=dict(genome_id=genome,distinct_proteins=g.protein.nunique(),same_contig_strand_triplet=bool(cand),local_triplet=any(c['within_declared_window'] for c in cand),local_triplet_all_partial00=any(c['within_declared_window'] and c['all_partial00'] for c in cand));out.append(rec)
 if cand:triples.append(cand[0])
a=pd.DataFrame(out);a.to_csv(O/'opu_organization_summary.tsv',sep='\t',index=False);pd.DataFrame(triples).to_csv(O/'opu_best_local_triplets.tsv',sep='\t',index=False)
res=[]
for q in ['primary','strict']:
 z=pd.read_csv(O/f'{q}_expanded_resource_calls.tsv',sep='\t').merge(a,on='genome_id',how='left',validate='one_to_one');z['local_triplet']=z.local_triplet.eq(True);z['local_triplet_all_partial00']=z.local_triplet_all_partial00.eq(True);z.to_csv(O/f'{q}_resource_locus_checked.tsv',sep='\t',index=False)
 for key in ['Opu_ABC_K05845_6_7','local_triplet','local_triplet_all_partial00']:
  res.append(dict(quality=q,definition=key,n=len(z),detected=int(z[key].sum()),samples=int(z.loc[z[key],'sample_id'].nunique()),regions=int(z.loc[z[key],'region_group'].nunique())))
pd.DataFrame(res).to_csv(O/'resource_local_triplet_summary.tsv',sep='\t',index=False);print(a.drop(columns='genome_id').sum().to_string());print(pd.DataFrame(res).to_string(index=False))
