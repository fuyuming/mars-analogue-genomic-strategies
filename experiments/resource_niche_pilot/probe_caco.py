"""Execute isolated upstream parsers on declared synthetic, non-biological fixtures."""
from pathlib import Path
import ast,os,json,hashlib
import pandas as pd
H=Path(__file__).resolve().parent;O=H/'host/parser_probe';O.mkdir(exist_ok=True)
# HMMER3 domtblout: target protein length field 2; query HMM length field 5;
# HMM start/end fields 15/16; alignment protein start/end fields 17/18.
def row(protein,tlen,family,qlen,hs,he,as_,ae,score=100):
 return ' '.join(map(str,[protein,'-',tlen,family,'-',qlen,'1e-30',score,0,1,1,'1e-30','1e-30',score,0,hs,he,as_,ae,as_,ae,.99,'synthetic']))+'\n'
fixtures={
 'long_protein_complete_domain':[row('p1',1000,'GH1.hmm',100,1,100,301,400)],
 'short_fragment_partial_domain':[row('p1',40,'GH1.hmm',400,1,40,1,40)],
 'two_families_one_protein':[row('p1',300,'GH1.hmm',150,1,150,1,150),row('p1',300,'CBM1.hmm',150,1,150,151,300,90)],
 'profile_order_overwrite':[row('p1',300,'GH1.hmm',150,1,150,1,150),row('p2',300,'GH1.hmm',150,1,150,1,150),row('p1',300,'CBM1.hmm',150,1,150,151,300,90)],
 'last_target_score_gate':[row('p1',100,'GH1.hmm',100,1,100,1,100,10)]}
results={}
for version,fn in [('archived','process_hmmsearch_output'),('current','parse_hmmsearch_output')]:
 src=H/'upstream'/version/'CaCo.py';tree=ast.parse(src.read_text());f=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==fn)
 ns={'pd':pd,'os':os};exec(compile(ast.Module(body=[f],type_ignores=[]),str(src),'exec'),ns)
 res={}
 for name,lines in fixtures.items():
  p=O/(version+'_'+name+'.domtbl');p.write_text('# synthetic fixture; NOT biological evidence\n'+''.join(lines));out=Path(str(p)+'.parsed')
  if out.exists():out.unlink()
  ns[fn](str(p),str(O));z=pd.read_csv(out,sep='\t');res[name]=z.to_dict('records')
 results[version]=res
assert results['archived']['long_protein_complete_domain']==[] and results['current']['long_protein_complete_domain']==[]
assert results['current']['short_fragment_partial_domain'][0]['Hit_Coverage']==1
assert len(results['current']['two_families_one_protein'])==1
assert next(x for x in results['current']['profile_order_overwrite'] if x['Query_Domain']=='p1')['Hit_Name']=='CBM1'
assert len(results['archived']['last_target_score_gate'])==1 and not results['current']['last_target_score_gate']
key=json.loads((H/'upstream/archived/data__substrate_key.json').read_text());subs={x.strip() for v in key.values() for x in v.split(',') if x.strip()}
classes={p:sum(k.startswith(p) for k in key) for p in ['AA','CBM','CE','GH','PL','GT']}
report={'synthetic_not_biological':True,'upstream_results':results,'model_coverage_expectations':{'long_protein_complete_domain':1.0,'short_fragment_partial_domain':.1},'HMMER_zero_based_field_semantics':{'target_protein_length':2,'query_model_length':5,'hmm_from':15,'hmm_to':16,'ali_from':17,'ali_to':18},'mapped_families':len(key),'unique_mapped_substrates':len(subs),'mapped_family_classes':classes,'mapping_identical_between_versions':(H/'upstream/archived/data__substrate_key.json').read_bytes()==(H/'upstream/current/data__substrate_key.json').read_bytes(),'limits':'Confirms parser behavior on explicit fixtures; not evidence that whole published ecological results are invalid. Biological impact assessed separately on the same pilot raw domtbl.'}
(O/'audit_checks.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items() if k!='upstream_results'},indent=2))
