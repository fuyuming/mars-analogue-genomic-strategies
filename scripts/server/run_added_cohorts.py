#!/usr/bin/env python3
"""Finite, checkpointed processing of a frozen supplementary MAG manifest.

Reuse existing versioned databases. Never restart or overwrite an unverified
analysis directory. This stops after quality/annotation/taxonomy, before ANI or
ecological admission. Host configuration is supplied separately, never embedded.
"""
import argparse, concurrent.futures, csv, fcntl, gzip, hashlib, json, os
from pathlib import Path
import subprocess, time

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def writejson(path,obj):
    p=Path(path);t=p.with_suffix(p.suffix+'.part');t.write_text(json.dumps(obj,indent=2));t.replace(p)

def table(path):
    with Path(path).open() as f:return list(csv.DictReader(f,delimiter='\t'))

def records(path):
    name=None;seq=[]
    with Path(path).open() as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].split()[0];seq=[]
            else:seq.append(line.strip())
    if name is not None:yield name,''.join(seq)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--config',type=Path,required=True)
    a=ap.parse_args();cfg=json.loads(a.config.read_text());root=Path(cfg['workdir'])
    root.mkdir(exist_ok=True,parents=True)
    lock=(root/'added_pipeline.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    lock.write(str(os.getpid()));lock.flush()
    manifest=Path(cfg['manifest']);rows=table(manifest);fingerprint=digest(manifest)
    assert len(rows)==cfg['expected_count'] and len({r['filename'] for r in rows})==len(rows)
    def state(stage,**kw):writejson(root/'pipeline_status.json',dict(stage=stage,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),pid=os.getpid(),manifest_sha256=fingerprint,**kw))
    def run(cmd,log,env=None):
        with Path(log).open('a') as f:
            f.write(repr(list(map(str,cmd)))+'\n');f.flush()
            subprocess.run(list(map(str,cmd)),stdout=f,stderr=f,check=True,env=env)
    def env_for(directory):
        e=os.environ.copy();e['PATH']=str(Path(directory)/'bin')+':'+e['PATH'];return e
    try:
        # Wait for the existing downloader's final receipt, without competing for files.
        state('waiting_for_existing_download')
        deadline=time.monotonic()+cfg.get('download_wait_hours',24)*3600
        oldstatus=root/'genomes/download_status.json'
        while True:
            try:s=json.loads(oldstatus.read_text())
            except (FileNotFoundError,json.JSONDecodeError):s={}
            if 'finished' in s:
                assert s['final_incomplete_or_bad']==0, 'Existing download failed; diagnose before any restart'
                break
            if time.monotonic()>deadline:raise TimeoutError('Existing downloader did not finish within finite wait')
            time.sleep(30)
        state('downloading_corrected_three')
        if not (root/'corrected_download_verified.json').exists():
            run(['python3',root/'download_mag_fastas.py','--manifest',root/'manifests/download_three.tsv','--outdir',root/'genomes','--parallel','3','--status-json',root/'genomes/corrected_three_status.json','--receipt',root/'genomes/corrected_three_receipt.tsv'],root/'corrected_three_wrapper.log')
            writejson(root/'corrected_download_verified.json',{'status':'finished'})
        inputs=root/'uniform_inputs';inputs.mkdir(exist_ok=True);aliases=[]
        state('verifying_all_inputs')
        for r in rows:
            p=root/'genomes'/r['filename'];assert p.stat().st_size==int(r['bytes'])
            data=p.read_bytes();assert hashlib.md5(data).hexdigest()==r['md5']
            raw=gzip.decompress(data);assert raw.startswith(b'>')
            gid='NEE56_'+r['filename'].removesuffix('.fna.gz')
            dest=inputs/(gid+'.fna')
            if dest.exists():assert digest(dest)==hashlib.sha256(raw).hexdigest()
            else:dest.write_bytes(raw)
            aliases.append(dict(genome_id=gid,compressed_file=r['filename'],source_url=r['url'],md5=r['md5'],sha256=hashlib.sha256(data).hexdigest(),uncompressed_sha256=digest(dest)))
        with (root/'input_aliases.tsv').open('w') as f:
            w=csv.DictWriter(f,aliases[0].keys(),delimiter='\t');w.writeheader();w.writerows(aliases)
        writejson(root/'input_verification.json',{'verified':len(aliases),'manifest_sha256':fingerprint,'compressed_bytes':sum(int(r['bytes']) for r in rows)})
        threads=int(cfg.get('threads',16));prod=Path(cfg['gtdb_env'])/'bin/prodigal'
        v=subprocess.run([prod,'-v'],capture_output=True,text=True,check=True);assert 'V2.6.3' in v.stdout+v.stderr
        state('prodigal_meta',genomes=len(aliases),workers=threads)
        def predict(r):
            gid=r['genome_id'];out=root/'prodigal'/gid;done=out/'verified.json'
            if done.exists():
                d=json.loads(done.read_text());assert d['input_sha256']==r['uncompressed_sha256']
                for n,h in d['outputs'].items():assert digest(out/n)==h
                return d
            assert not out.exists(),f'Unverified prediction directory {gid}; inspect before rerun'
            out.mkdir(parents=True)
            cmd=[prod,'-i',inputs/(gid+'.fna'),'-a',out/'proteins.faa','-d',out/'genes.ffn','-o',out/'genes.gff','-f','gff','-p','meta','-q']
            run(cmd,out/'run.log');aa=list(records(out/'proteins.faa'));nt=list(records(out/'genes.ffn'))
            cds=sum('\tCDS\t' in line for line in (out/'genes.gff').open())
            assert len(aa)==len(nt)==cds and cds>0 and len({n for n,s in aa})==cds
            assert 'run_type=Metagenomic' in (out/'genes.gff').read_text()[:1000]
            kept=0
            with (out/'prefixed_min20.faa').open('w') as f:
                for n,s in aa:
                    s=s.replace('*','')
                    if len(s)>=20:f.write('>'+gid+'|'+n+'\n'+s+'\n');kept+=1
            d=dict(genome_id=gid,input_sha256=r['uncompressed_sha256'],proteins=cds,retained_min20=kept,prodigal_binary_sha256=digest(prod),outputs={n:digest(out/n) for n in ['proteins.faa','genes.ffn','genes.gff','prefixed_min20.faa']})
            writejson(done,d);return d
        with concurrent.futures.ThreadPoolExecutor(threads) as ex:pred=list(ex.map(predict,aliases))
        writejson(root/'prediction_summary.json',{'verified':len(pred),'proteins':sum(x['proteins'] for x in pred),'retained_min20':sum(x['retained_min20'] for x in pred)})
        qc=root/'checkm2_1.1.0';qc_done=root/'checkm2_verified.json';state('checkm2',threads=threads)
        if not qc_done.exists():
            assert not qc.exists(),'Unverified CheckM2 directory; inspect before rerun'
            run([Path(cfg['checkm_env'])/'bin/checkm2','predict','--input',inputs,'--output-directory',qc,'--threads',threads,'--extension','fna','--database_path',cfg['checkm_db']],root/'checkm2_wrapper.log',env_for(cfg['checkm_env']))
            q=table(qc/'quality_report.tsv');assert {r['Name'] for r in q}=={r['genome_id'] for r in aliases} and len(q)==len(aliases)
            writejson(qc_done,{'sha256':digest(qc/'quality_report.tsv'),'n':len(q)})
        else:assert digest(qc/'quality_report.tsv')==json.loads(qc_done.read_text())['sha256']
        q=table(qc/'quality_report.tsv');accepted=[];strict=[]
        for r in q:
            r['primary']=float(r['Completeness'])>=50 and float(r['Contamination'])<10
            r['strict']=float(r['Completeness'])>=90 and float(r['Contamination'])<5
            if r['primary']:accepted.append(r['Name'])
            if r['strict']:strict.append(r['Name'])
        assert accepted
        with (root/'cohort_quality_registry.tsv').open('w') as f:
            w=csv.DictWriter(f,q[0].keys(),delimiter='\t');w.writeheader();w.writerows(q)
        writejson(root/'quality_gate_summary.json',dict(total=len(q),primary=len(accepted),strict=len(strict),thresholds='completeness>=50,contamination<10; strict>=90,<5',ecological_admission=False))
        primary=root/'primary_inputs';primary.mkdir(exist_ok=True)
        for gid in accepted:
            dst=primary/(gid+'.fna')
            if not dst.exists():dst.symlink_to(inputs/(gid+'.fna'))
        combined=root/'primary_uniform.faa'
        if not combined.exists():
            with combined.with_suffix('.part').open('wb') as f:
                for gid in sorted(accepted):f.write((root/'prodigal'/gid/'prefixed_min20.faa').read_bytes())
            combined.with_suffix('.part').replace(combined)
        K=Path(cfg['kofam_dir']);G=Path(cfg['gtdb_data'])
        writejson(root/'annotation_provenance.json',dict(combined_sha256=digest(combined),kofam_config_sha256=digest(K/'config.yml'),ko_list_sha256=digest(K.parent/'ko_list'),prokaryote_hal_sha256=digest(K.parent/'profiles/prokaryote.hal'),gtdb_metadata_sha256=digest(G/'metadata/metadata.txt'),checkm_db_bytes=Path(cfg['checkm_db']).stat().st_size))
        state('kofam',primary=len(accepted),strict=len(strict),threads=threads)
        ko=root/'full_kofam_1.3.0';kd=root/'kofam_verified.json'
        if not kd.exists():
            assert not ko.exists(),'Unverified Kofam directory; inspect before rerun'
            ko.mkdir();result=ko/'mapper_one_line.tsv'
            run([K/'exec_annotation','-c',K/'config.yml','--cpu',threads,'--tmp-dir',ko/'tmp','--format','mapper-one-line','--no-report-unannotated','-o',result,combined],ko/'run.log')
            allids={n for n,s in records(combined)};ids=set()
            for line in result.open():
                if line.strip() and not line.startswith('#'):ids.add(line.split('\t')[0])
            assert ids<=allids
            writejson(kd,dict(annotated=len(ids),input=len(allids),sha256=digest(result)))
        state('gtdb_classify',threads=threads,primary=len(accepted))
        gout=root/'gtdbtk_r226';gd=root/'gtdb_verified.json'
        if not gd.exists():
            assert not gout.exists(),'Unverified GTDB directory; inspect before rerun'
            e=env_for(cfg['gtdb_env']);e['GTDBTK_DATA_PATH']=str(G)
            run([Path(cfg['gtdb_env'])/'bin/gtdbtk','classify_wf','--genome_dir',primary,'--out_dir',gout,'--extension','fna','--cpus',threads,'--pplacer_cpus',threads],root/'gtdb_wrapper.log',e)
            ids=[]
            for name in ['gtdbtk.bac120.summary.tsv','gtdbtk.ar53.summary.tsv']:
                p=gout/'classify'/name
                if p.exists():ids += [r['user_genome'] for r in table(p)]
            assert set(ids)==set(accepted) and len(ids)==len(accepted),'Taxonomy not complete; inspect failed-genome outputs'
            writejson(gd,dict(classified=len(ids),expected=len(accepted)))
        state('quality_annotation_taxonomy_complete',primary=len(accepted),strict=len(strict),note='ANI, source/location deduplication, ecological admission and expanded tree pending')
    except Exception as ex:state('failed',error=repr(ex));raise

if __name__=='__main__':main()
