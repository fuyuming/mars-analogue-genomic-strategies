"""Reuse frozen markers and biosynthesis counts; no gene reannotation."""
from pathlib import Path
import os,json,hashlib
import pandas as pd
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=Path(os.environ.get('NEE_PROJECT_ROOT',str(HERE.parents[1])))
RAW=ROOT.parent/'code/result_raw'; PREV=ROOT/'NEE_revision_20260926/14_biosynthetic_strategy_20260927/results'
OUT=HERE/'inputs';OUT.mkdir(exist_ok=True)
traits=['biofilm_eps_surface','osmotic_desiccation_salt','oxidative_redox','trace_gas_energy',
 'dormancy_resuscitation','dna_repair_radiation','cold_protein_quality','sulfur_chemolithotrophy','AA_count','cofactor_count']
denom=[4,15,9,9,4,12,3,7,12,8]
paths=[RAW/'KO_panel_repair_20260926/amended_inputs/core65_MAG_module_breadth.tsv.gz',PREV/'genome_route_traits.tsv',
 PREV/'Qaidam_sample_traits.tsv',ROOT/'NEE_revision_20260926/07_scene_function_analysis_20260927/sample_and_spore_inputs/qaidam_MAG_to_sample_core65.tsv']
module=pd.read_csv(paths[0],sep='\t');genome=pd.read_csv(paths[1],sep='\t').set_index('genome_id')
counts=module.pivot(index='genome_id',columns='module',values='n_detected_KOs').astype(int)
d=genome.join(counts,validate='one_to_one')
assert len(d)==3904
b=d[d.gtdb_domain.eq('Bacteria') & d.species_representative & d.passes_strict_90_5].copy()
assert len(b)==879 and b.light_energy.eq(0).all()
cov=['checkm2_completeness','checkm2_contamination','log2_genome_size']
b['log2_genome_size']=np.log2(b.genome_size)
means=b[cov].mean();sds=b[cov].std(ddof=0)
for c in cov:b['z_'+c]=(b[c]-means[c])/sds[c]
keep=['dataset_id','family_path','ani_cluster_95']+cov+['z_'+c for c in cov]+traits
b[keep].reset_index().to_csv(OUT/'global_strict.tsv',sep='\t',index=False)
shared=b.groupby('family_path').dataset_id.nunique();b[b.family_path.isin(shared[shared>=2].index)][keep].reset_index().to_csv(OUT/'global_shared_families.tsv',sep='\t',index=False)

mapping=pd.read_csv(paths[3],sep='\t')[['genome_id','sample_id','site_id']]
q=d[d.gtdb_domain.eq('Bacteria')].reset_index().merge(mapping,on='genome_id',validate='one_to_one');assert len(q)==1492
samples=pd.read_csv(paths[2],sep='\t')
env_scales={}
for name,oldname,strict in [('primary','primary_min5',False),('strict','strict90_5_min3',True)]:
 s=samples[samples.stratum.eq(oldname)].copy().reset_index(drop=True)
 mag=q[q.passes_strict_90_5] if strict else q
 k=mag.groupby('sample_id')[traits].sum().loc[s.sample_id]
 nmag=mag.groupby('sample_id').size().loc[s.sample_id]
 assert np.array_equal(nmag.to_numpy(),s.n_MAGs.to_numpy())
 for trait in traits:s[trait]=k[trait].to_numpy()
 s['log2_genome_size']=np.log2(s.genome_size)
 envcols=[]
 for exposure in ['water','EC','TOC']:
  base='log2_'+exposure
  s[base+'_between']=s.groupby('site_id')[base].transform('mean')
  s[base+'_within']=s[base]-s[base+'_between']
  envcols.extend([base+'_within',base+'_between'])
 controls=['pH','depth_cm','checkm2_completeness','checkm2_contamination','log2_genome_size','family_PC1','family_PC2','family_PC3']
 mm=s[envcols+controls].mean();ss=s[envcols+controls].std(ddof=0)
 assert (ss>1e-8).all()
 for c in envcols+controls:s['z_'+c]=(s[c]-mm[c])/ss[c]
 s.to_csv(OUT/f'qaidam_{name}.tsv',sep='\t',index=False)
 env_scales[name]={'means':mm.to_dict(),'sds':ss.to_dict(),'environment_columns':envcols,'control_columns':controls,
  'n_samples':len(s),'n_sites':s.site_id.nunique(),'n_MAGs':int(s.n_MAGs.sum()),
  'within_varying_sites':{e:int((s.groupby('site_id')['log2_'+e].nunique()>1).sum()) for e in ['water','EC','TOC']}}
meta={'traits':traits,'denominators':denom,'global_covariates':cov,'global_means':means.to_dict(),'global_sds':sds.to_dict(),
 'environment_scales':env_scales,'global_n':len(b),'source_counts':b.groupby('dataset_id').size().to_dict(),
 'inputs':[{'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]}
(OUT/'metadata.json').write_text(json.dumps(meta,indent=2))
print(json.dumps({k:meta[k] for k in ['global_n','source_counts','environment_scales']},indent=2))
