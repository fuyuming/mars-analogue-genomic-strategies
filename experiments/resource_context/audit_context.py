from pathlib import Path
import csv,json,hashlib,collections
H=Path(__file__).resolve().parent;prev=H/'inputs';(H/'outputs').mkdir(exist_ok=True)
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def put(name,rows):
 with (H/'outputs'/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
ctx=tab(prev/'candidate_context_verified_KOs.tsv');receipt=json.loads((H/'inputs/context_KO_complete.json').read_text());assert hashlib.sha256((prev/'candidate_context_verified_KOs.tsv').read_bytes()).hexdigest()==receipt['output_sha256']
ko={r['knum']:r['definition'] for r in tab(H/'inputs/ko_definitions.tsv')};cand=tab(prev/'conservative_substrate_candidates.tsv');lookup={r['protein']:r for r in cand};assert len(lookup)==57
out=[];detailed=[]
for target,c in lookup.items():
 window=[r for r in ctx if r['target']==target];selfrow=next(r for r in window if int(r['offset'])==0);neighbors=[r for r in window if int(r['offset'])!=0]
 annotated=[];transport=[]
 for r in window:
  defs='; '.join(k+': '+ko[k] for k in r['accepted_KOs'].split(';') if k);dist=max(0,int(r['start'])-int(selfrow['end'])-1,int(selfrow['start'])-int(r['end'])-1)
  # Lexical flags locate records for manual review; do not assign uptake specificity.
  flag=bool(defs) and any(x in defs.lower() for x in ['transport','permease','symporter','antiporter','porin'])
  detailed.append(dict(**r,KO_definitions=defs,same_strand=r['strand']==selfrow['strand'],intergenic_distance_to_target_bp=dist,transport_text_flag=flag))
  if r['offset']!='0' and defs:annotated.append(r['offset']+' '+defs)
  if r['offset']!='0' and flag:transport.append(r['offset']+' '+defs)
 out.append(dict(target=target,genome=c['genome'],family=c['anchor_family'],substrate_label=c['labels'],target_KOs=selfrow['accepted_KOs'],target_definitions='; '.join(ko[k] for k in selfrow['accepted_KOs'].split(';') if k),actual_neighbors=len(neighbors),neighbors_with_KO=sum(bool(r['accepted_KOs']) for r in neighbors),near_ORF_end=selfrow['target_near_contig_ORF_end'],transport_text_flag_neighbors=' || '.join(transport),annotated_neighbors=' || '.join(annotated)))
put('candidate_context_annotated.tsv',detailed);put('candidate_target_audit.tsv',out)
windowp={r['protein'] for r in ctx};targets=set(lookup);non_targets=windowp-targets
summary=dict(targets=len(targets),window_records=len(ctx),unique_window_proteins=len(windowp),targets_with_KO=sum(bool(r['target_KOs']) for r in out),unique_non_target_proteins=len(non_targets),non_target_proteins_with_KO=len({r['protein'] for r in ctx if r['protein'] in non_targets and r['accepted_KOs']}),targets_with_at_least_one_annotated_neighbor=sum(r['neighbors_with_KO']>0 for r in out),targets_with_transport_text_flag=sum(bool(r['transport_text_flag_neighbors']) for r in out),target_near_ORF_end=sum(r['near_ORF_end']=='True' for r in out),labels=dict(collections.Counter(r['substrate_label'] for r in out)),note='Window totals include target proteins. Text flags require biological review; no uptake specificity inferred.')
(H/'outputs/context_audit_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
for r in out:
 print(r['target'].split('|')[-1],r['family'],r['substrate_label'],'SELF:',r['target_definitions'],'TRANSPORT:',r['transport_text_flag_neighbors'])
