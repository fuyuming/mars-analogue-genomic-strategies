from pathlib import Path
import csv,json,collections,hashlib
H=Path(__file__).resolve().parent;REV=H.parent;(H/'outputs').mkdir(exist_ok=True)
RULES={'KdpABC':['K01546','K01547','K01548'],'ProVWX':['K02000','K02001','K02002'],'EctABC':['K06718','K00836','K06720'],'BetAB':['K00108','K00130'],'CoxCut':['K03518','K03519','K03520']}
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def put(n,rs):
 with (H/'outputs'/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
receipt=json.loads((H/'inputs/export_receipt.json').read_text());assert hashlib.sha256((H/'inputs/new556_restored74_presence.tsv').read_bytes()).hexdigest()==receipt['output_sha256']
mat={r['genome_id']:r for r in tab(H/'inputs/new556_restored74_presence.tsv')};manifest=tab(REV/'final_phylogeny/frozen_genome_manifest.tsv');states=[];summaries=[];gates=[];familycells=[];linksummary=[]
for quality in ['primary','strict']:
 links={r['Name']:r for r in tab(H/f'inputs/source_checked_{quality}_resource_links.tsv')}
 selected=[r for r in manifest if r[quality+'_representative']=='True' and r['genome_id'] in mat]
 for r in selected:
  row=mat[r['genome_id']];family=next((x for x in r['taxonomy_r226'].split(';') if x.startswith('f__')),'f__')
  for name,kos in RULES.items():
   n=sum(int(row[k]) for k in kos);state='all_detected' if n==len(kos) else 'partial' if n else 'not_detected'
   states.append(dict(quality=quality,genome_id=r['genome_id'],dataset_id=r['dataset_id'],family=family,origin_unit=r['origin_unit'],role=r['role'],resource_eligible=r['genome_id'] in links,system=name,state=state,detected=n,required=len(kos)))
 for dataset in sorted({r['dataset_id'] for r in selected}):
  ss=[r for r in selected if r['dataset_id']==dataset];families={r['genome_id']:next((x for x in r['taxonomy_r226'].split(';') if x.startswith('f__')),'f__') for r in ss}
  for name,kos in RULES.items():
   c=collections.Counter(next(s['state'] for s in states if s['quality']==quality and s['genome_id']==r['genome_id'] and s['system']==name) for r in ss)
   summaries.append(dict(quality=quality,dataset=dataset,system=name,n=len(ss),all_detected=c['all_detected'],partial=c['partial'],not_detected=c['not_detected']))
  for maint in list(RULES)[:-1]:
   eligible={g:[] for g in [2,3,5]};state_supported={g:[] for g in [2,3,5]}
   for fam in sorted(set(families.values())-{'f__'}):
    members=[r for r in ss if families[r['genome_id']]==fam];counts=collections.Counter();units=collections.Counter()
    for r in members:
     row=mat[r['genome_id']];m=all(int(row[k]) for k in RULES[maint]);e=all(int(row[k]) for k in RULES['CoxCut']);counts[2*int(m)+int(e)]+=1
     # Only unambiguous single-parent, source-checked samples count as independent units.
     if r['genome_id'] in links:units[links[r['genome_id']]['sample_id']]+=1
    cells=[counts[j] for j in range(4)];familycells.append(dict(quality=quality,dataset=dataset,pair=maint+'__CoxCut',family=fam,n=len(members),n00=cells[0],n01=cells[1],n10=cells[2],n11=cells[3],source_checked_samples_ge3=sum(v>=3 for v in units.values())))
    for g in eligible:
     if min(cells)>=g:state_supported[g].append(fam)
     if min(cells)>=g and sum(v>=3 for v in units.values())>=2:eligible[g].append(fam)
   for g,fams in eligible.items():gates.append(dict(quality=quality,dataset=dataset,pair=maint+'__CoxCut',min_per_cell=g,families_with_all_four_states_ge_gate=len(state_supported[g]),state_supported_families=';'.join(state_supported[g]),source_checked_Mackay_eligible_families=len(fams) if dataset=='mackay_PRJNA630822' else 'not_evaluated',families=';'.join(fams)))
 mm=[r for r in selected if r['genome_id'] in links];linksummary.append(dict(quality=quality,all_source_checked_MAGs=len(links),global_representatives_with_checked_resource=len(mm),samples=len({links[r['genome_id']]['sample_id'] for r in mm}),regions=len({links[r['genome_id']]['region_group'] for r in mm})))
put('new_representative_system_states.tsv',states);put('new_representative_system_summary.tsv',summaries);put('joint_support_feasibility.tsv',gates);put('family_joint_cells.tsv',familycells);put('resource_units_after_global_representative_selection.tsv',linksummary)
print(json.dumps(linksummary,indent=2));print('MODULES');print('\n'.join(str(r) for r in summaries));print('GATE3');print('\n'.join(str(r) for r in gates if r['min_per_cell']==3))
