from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
H=Path(__file__).resolve().parent;R=H.parent;O=H/'host';O.mkdir(exist_ok=True)
paths=[O/'frozen/all5477_restored74_presence.tsv.gz',O/'frozen/union_genome_manifest.tsv',O/'frozen/final_tree_tip_manifest.tsv']
x=pd.read_csv(paths[0],sep='\t',index_col=0);meta=pd.read_csv(paths[1],sep='\t').set_index('tip_id');tips=pd.read_csv(paths[2],sep='\t')
for rank,pref in [('phylum','p__'),('class','c__'),('order','o__'),('family','f__'),('genus','g__')]:meta[rank]=meta.taxonomy_r226.map(lambda t:next((z[len(pref):] for z in t.split(';') if z.startswith(pref)),'') or 'Unclassified')
signatures={'KdpABC':['K01546','K01547','K01548'],'KdpDE':['K07646','K07667'],'ProVWX':['K02000','K02001','K02002'],'BetAB':['K00108','K00130']};pairs=[('KdpABC','ProVWX'),('BetAB','ProVWX')];rows=[];counts=[];fc=[];pc=[];resources=[];resourcecounts=[];anomalies=[]
for quality in ['primary','strict']:
 m=meta.loc[tips.loc[tips.quality.eq(quality),'tip_id']].copy();a=x.loc[m.genome_id].copy();a.index=m.index;assert a.notna().all().all();m['quality']=quality
 for name,kos in signatures.items():
  m[name+'_n']=a[kos].sum(axis=1);m[name+'_status']=np.select([m[name+'_n'].eq(0),m[name+'_n'].eq(len(kos))],['none','all_listed'],'partial');m[name]=m[name+'_n'].eq(len(kos))
 for domain,d in m.groupby('expected_domain'):
  for name,kos in signatures.items():
   c=d[name+'_status'].value_counts();counts.append(dict(quality=quality,domain=domain,signature=name,denominator=len(d),none=int(c.get('none',0)),partial=int(c.get('partial',0)),all_listed=int(c.get('all_listed',0))))
  for A,B in pairs:
   st=np.select([d[A]&d[B],d[A]&~d[B],~d[A]&d[B]],['both','A_only','B_only'],'neither');d=d.copy();d['pair_state']=st
   for state,n in d.pair_state.value_counts().items():pc.append(dict(quality=quality,domain=domain,pair=A+'__'+B,state=state,n=n))
   for fam,g in d[d.family.ne('Unclassified')].groupby('family'):
    c=g.pair_state.value_counts();ca=g[g.pair_state.eq('A_only')].dataset_id.nunique();cb=g[g.pair_state.eq('B_only')].dataset_id.nunique();na=int(c.get('A_only',0));nb=int(c.get('B_only',0));cross=sum(gg.pair_state.isin(['A_only','B_only']).sum()>0 and {'A_only','B_only'}.issubset(set(gg.pair_state)) for _,gg in g.groupby('dataset_id'))
    fc.append(dict(quality=quality,domain=domain,pair=A+'__'+B,family=fam,n=len(g),A_only=na,B_only=nb,both=int(c.get('both',0)),neither=int(c.get('neither',0)),A_only_catalogues=ca,B_only_catalogues=cb,within_catalogue_both_states=cross,coarse_comparability=(na>=3 and nb>=3 and ca>=2 and cb>=2)))
  for flag,idx in [('KdpDE_all_without_KdpABC_all',d.KdpDE&~d.KdpABC),('KdpABC_all_without_KdpDE_all',d.KdpABC&~d.KdpDE)]:anomalies.append(dict(quality=quality,domain=domain,configuration=flag,n=int(idx.sum())))
 rows.append(m.reset_index().drop(columns=['fasta_path'],errors='ignore'))
 rp=O/f'frozen/source_checked_{quality}_resource_links.tsv';paths.append(rp);rl=pd.read_csv(rp,sep='\t');z=m.reset_index().merge(rl[['genome_id','sample_id','region_group','water_percent','EC_dS_m','pH','TOC_value','TOC_censored']],on='genome_id',how='inner',validate='one_to_one');resources.append(z.drop(columns=['fasta_path'],errors='ignore'))
 for name in signatures:
  resourcecounts.append(dict(quality=quality,signature=name,MAGs=len(z),samples=z.sample_id.nunique(),regions=z.region_group.nunique(),all_listed=int(z[name].sum()),samples_with_all=z.loc[z[name],'sample_id'].nunique(),samples_with_not_all=z.loc[~z[name],'sample_id'].nunique(),families_with_both_states=int(z.groupby('family')[name].nunique().eq(2).sum())))
pd.concat(rows).to_csv(O/'signature_calls.tsv',sep='\t',index=False);pd.DataFrame(counts).to_csv(O/'signature_counts.tsv',sep='\t',index=False);pd.DataFrame(pc).to_csv(O/'pair_states.tsv',sep='\t',index=False);f=pd.DataFrame(fc);f.to_csv(O/'family_comparability.tsv',sep='\t',index=False);f[f.coarse_comparability].to_csv(O/'comparable_families.tsv',sep='\t',index=False);pd.DataFrame(anomalies).to_csv(O/'pump_regulator_discordance.tsv',sep='\t',index=False);pd.concat(resources).to_csv(O/'source_checked_resource_signature_calls.tsv',sep='\t',index=False);pd.DataFrame(resourcecounts).to_csv(O/'resource_identifiability.tsv',sep='\t',index=False)
(O/'input_manifest.json').write_text(json.dumps([dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],indent=2));print(pd.DataFrame(counts).to_string(index=False));print('COMPARABLE',f[f.coarse_comparability].to_string(index=False));print(pd.DataFrame(resourcecounts).to_string(index=False))
