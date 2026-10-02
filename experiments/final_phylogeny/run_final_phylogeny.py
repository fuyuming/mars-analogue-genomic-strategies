"""Uniform GTDB markers and independently inferred quality-specific phylogenies.
No overwrite/restart of unverified stages. Existing completed outputs are hash checked.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
from collections import Counter
import argparse,csv,json,hashlib,subprocess,time,gzip,fcntl,os,traceback,re

def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for chunk in iter(lambda:f.read(4*1024*1024),b''):h.update(chunk)
 return h.hexdigest()
def write_json(p,d):
 p=Path(p);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(d,indent=2));tmp.replace(p)
def table(p):return list(csv.DictReader(Path(p).open(),delimiter='\t'))
def write_table(p,rows,fields=None):
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def read_fasta(p):
 out={};key=None;op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt') as f:
  for line in f:
   if line.startswith('>'):
    key=line[1:].split()[0];assert key not in out,('duplicate ID',key);out[key]=''
   elif line.strip():
    assert key is not None;out[key]+=line.strip().upper()
 return out

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True);args=ap.parse_args();cfg=json.loads(args.config.read_text())
 R=Path(cfg['root']);R.mkdir(exist_ok=True,parents=True);lock=(R/'pipeline.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 manifest=R/'union_genome_manifest.tsv';rows=table(manifest);expected={r['tip_id']:r for r in rows};fingerprint=digest(manifest);assert len(expected)==cfg['expected_union'] and len({r['genome_id'] for r in rows})==len(rows)
 assert fingerprint==cfg['manifest_sha256'];E=Path(cfg['gtdb_env']);DB=Path(cfg['gtdb_data']);env=os.environ.copy();env['PATH']=str(E/'bin')+':'+env['PATH'];env['GTDBTK_DATA_PATH']=str(DB)
 def state(stage,**kw):
  data=dict(stage=stage,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),manifest_sha256=fingerprint,**kw);write_json(R/'status.json',data);print(json.dumps(data),flush=True)
 def run(stage,cmd,artifacts,threads):
  done=R/(stage+'.done.json');cmd=list(map(str,cmd))
  if done.exists():
   d=json.loads(done.read_text());assert d['args']==cmd and d['manifest_sha256']==fingerprint
   for f,h in d['output_sha256'].items():assert digest(f)==h
   return
  assert not (R/(stage+'.started.json')).exists(),f'{stage} already attempted; inspect before explicit checkpoint resume'
  write_json(R/(stage+'.started.json'),dict(args=cmd,manifest_sha256=fingerprint,threads=threads,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
  t=time.monotonic()
  with (R/(stage+'.log')).open('x') as log:
   p=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,env=env,cwd=R)
   write_json(R/(stage+'.process.json'),dict(pid=p.pid,args=cmd));rc=p.wait()
  if rc!=0:raise RuntimeError((stage,'exit',rc))
  for p in artifacts:assert Path(p).is_file() and Path(p).stat().st_size>0,(stage,'missing output',str(p))
  write_json(done,dict(args=cmd,manifest_sha256=fingerprint,elapsed_seconds=time.monotonic()-t,exit=rc,output_sha256={str(p):digest(p) for p in artifacts}))
 try:
  state('input_preflight',expected_union=len(rows));genomes=R/'genomes';genomes.mkdir(exist_ok=True)
  def check(r):
   src=Path(r['fasta_path']);assert src.is_file() and src.stat().st_size>0,src
   with src.open('rb') as f:assert f.read(1)==b'>',('not plain FASTA',str(src))
   return dict(tip_id=r['tip_id'],genome_id=r['genome_id'],path=str(src),bytes=src.stat().st_size,sha256=digest(src))
  checks=[]
  with ThreadPoolExecutor(max_workers=8) as ex:
   for result in ex.map(check,rows):checks.append(result)
  for r in rows:
   target=genomes/(r['tip_id']+'.fna');src=Path(r['fasta_path'])
   if target.is_symlink():assert target.resolve()==src.resolve()
   else:assert not target.exists();target.symlink_to(src)
  assert len(list(genomes.glob('*.fna')))==len(rows)
  write_table(R/'input_fasta_checksums.tsv',checks)
  versions={}
  for name,cmd in [('gtdbtk',[E/'bin/gtdbtk','--version']),('iqtree',[cfg['iqtree'],'--version'])]:
   p=subprocess.run(list(map(str,cmd)),env=env,capture_output=True,text=True,check=True);versions[name]=p.stdout+p.stderr
  assert '2.5.2' in versions['gtdbtk'] and '2.0.7' in versions['iqtree']
  write_json(R/'input_preflight.json',dict(manifest_sha256=fingerprint,verified=len(checks),input_bytes=sum(x['bytes'] for x in checks),software=versions,gtdb_metadata_sha256=digest(DB/'metadata/metadata.txt'),resource_limits=dict(marker_threads=cfg['marker_threads'],concurrent_trees=2,tree_threads_each=cfg['tree_threads'],memory_each=cfg['tree_memory'])))
  state('identify',genomes=len(rows),threads=cfg['marker_threads'])
  run('identify',[E/'bin/gtdbtk','identify','--genome_dir',genomes,'--out_dir',R/'markers','--extension','fna','--cpus',cfg['marker_threads']],[],cfg['marker_threads'])
  state('align',threads=cfg['marker_threads'])
  paths={m:R/f'alignment/align/gtdbtk.{m}.user_msa.fasta.gz' for m in ['bac120','ar53']}
  run('align',[E/'bin/gtdbtk','align','--identify_dir',R/'markers','--out_dir',R/'alignment','--skip_gtdb_refs','--min_perc_aa','10','--cpus',cfg['marker_threads']],list(paths.values()),cfg['marker_threads'])
  state('alignment_QC')
  seqs={m:read_fasta(p) for m,p in paths.items()};assigned={};occupancy=[];issues=[]
  for marker,domain,L in [('bac120','Bacteria',5036),('ar53','Archaea',8062)]:
   for tip,s in seqs[marker].items():
    assert tip in expected and tip not in assigned,('unexpected/duplicate domain tip',tip)
    if len(s)!=L:issues.append(dict(tip_id=tip,issue='alignment_length',observed=len(s),expected=L))
    if set(s)-set('ACDEFGHIKLMNPQRSTVWYXBZJUO?-.'):issues.append(dict(tip_id=tip,issue='unexpected_symbols'))
    frac=sum(c in 'ACDEFGHIKLMNPQRSTVWY' for c in s)/L
    if frac<.10:issues.append(dict(tip_id=tip,issue='below_10_percent_effective_aa',observed=frac))
    want=expected[tip]['expected_domain']
    if want not in [domain,'Unclassified']:issues.append(dict(tip_id=tip,issue='domain_conflict',observed=domain,expected=want))
    assigned[tip]=marker;occupancy.append(dict(tip_id=tip,genome_id=expected[tip]['genome_id'],marker=marker,domain=domain,length=len(s),observed_aa_fraction=frac,original_domain=want))
  missing=set(expected)-set(assigned);filtered={}
  for marker,domain in [('bac120','Bacteria'),('ar53','Archaea')]:
   p=R/f'alignment/align/gtdbtk.{marker}.filtered.tsv'
   if not p.exists():continue
   for line in p.open():
    if not line.strip():continue
    fields=line.rstrip('\n').split('\t');tip=fields[0]
    reason=fields[1] if len(fields)==2 else ''
    match=re.fullmatch(r'Insufficient number of amino acids in MSA \((\d+(?:\.\d+)?)%\)',reason)
    # GTDB reports percentages rounded to one decimal; 10.0 can mean <10.
    if tip not in missing or tip in filtered or not match or not 0<=float(match.group(1))<=10.0:
     issues.append(dict(tip_id=tip,issue='unexpected_alignment_filter_record',file=str(p),reason=reason));continue
    if expected[tip]['expected_domain'] not in [domain,'Unclassified']:
     issues.append(dict(tip_id=tip,issue='filtered_domain_conflict',observed=domain,expected=expected[tip]['expected_domain']));continue
    filtered[tip]=str(p.relative_to(R))+': '+reason
  exclusions=[dict(tip_id=k,genome_id=expected[k]['genome_id'],reason=filtered[k]) for k in sorted(filtered)]
  for k in sorted(missing-set(filtered)):issues.append(dict(tip_id=k,issue='missing_without_documented_alignment_filter'))
  write_json(R/'alignment_QC.json',dict(manifest_sha256=fingerprint,input=len(rows),aligned=len(assigned),filtered_exclusions=exclusions,issues=issues,domains={m:len(d) for m,d in seqs.items()}));write_table(R/'marker_occupancy.tsv',occupancy);write_table(R/'marker_filter_exclusions.tsv',exclusions,['tip_id','genome_id','reason'])
  assert not issues,('Alignment QC failed; inspect alignment_QC.json',len(issues))
  tasks=[];final_rows=[];counts={}
  for quality in ['primary','strict']:
   counts[quality]={};out=R/'inference'/quality;out.mkdir(parents=True,exist_ok=True)
   for marker in ['bac120','ar53']:
    ids=sorted(k for k in seqs[marker] if expected[k][quality+'_representative']=='True');assert len(ids)>=4
    p=out/(marker+'.fasta');p.write_text(''.join('>'+k+'\n'+seqs[marker][k]+'\n' for k in ids))
    counts[quality][marker]=len(ids)
    duplicates=Counter(seqs[marker][k] for k in ids)
    write_json(out/(marker+'.alignment.json'),dict(tips=len(ids),sha256=digest(p),identical_alignment_groups=sum(v>1 for v in duplicates.values()),max_identical_group=max(duplicates.values())))
    for k in ids:final_rows.append(dict(quality=quality,marker=marker,tip_id=k,genome_id=expected[k]['genome_id'],ani_cluster_95=expected[k][quality+'_ANI95']))
    tasks.append((quality,marker,p,ids))
  write_table(R/'final_tree_tip_manifest.tsv',final_rows);write_json(R/'final_tree_counts.json',counts)
  state('tree_inference',counts=counts,max_total_threads=2*cfg['tree_threads'],memory_each=cfg['tree_memory'])
  def infer_quality(quality):
   for q,marker,msa,ids in tasks:
    if q!=quality:continue
    prefix=msa.parent/(marker+'_LG_F_R10');stage='iqtree_'+quality+'_'+marker
    outputs=[Path(str(prefix)+suffix) for suffix in ['.treefile','.iqtree','.contree','.log']]
    cmd=[cfg['iqtree'],'-s',msa,'-m','LG+F+R10','-B','1000','--alrt','1000','-T',cfg['tree_threads'],'--mem',cfg['tree_memory'],'-seed',cfg['seed'],'--prefix',prefix]
    run(stage,cmd,outputs,cfg['tree_threads'])
    found=re.findall(r'NEE60_G_\d{5}',outputs[0].read_text());assert len(found)==len(ids) and set(found)==set(ids),(stage,'tip mismatch')
    report=outputs[1].read_text();assert 'SH-aLRT' in report and '1000' in report and 'LG+F+R10' in report,(stage,'incomplete model/support report')
    write_json(R/(stage+'.validated.json'),dict(exact_tips=True,tips=len(ids),tree_sha256=digest(outputs[0]),report_sha256=digest(outputs[1]),note='Inference complete; biological branch/long-branch/old-tree comparison audit still required'))
  def guarded_infer(quality):
   write_json(R/(quality+'_tree_status.json'),dict(stage='running'))
   try:
    infer_quality(quality)
    write_json(R/(quality+'_tree_status.json'),dict(stage='complete'))
   except Exception as e:
    write_json(R/(quality+'_tree_status.json'),dict(stage='failed',error=repr(e)))
    raise
  with ThreadPoolExecutor(max_workers=2) as ex:
   futures=[ex.submit(guarded_infer,q) for q in ['primary','strict']]
   for f in as_completed(futures):f.result()
  state('complete',counts=counts,note='Independent quality-specific inference complete; final biological/topology/figure review pending')
 except Exception as e:
  state('failed',error=repr(e));traceback.print_exc();raise
if __name__=='__main__':main()
