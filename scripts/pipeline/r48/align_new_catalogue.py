"""Reuse finished annotations to align a new catalogue to corrected endpoints; no ecology fit."""
from pathlib import Path
import pandas as pd,json,hashlib
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/amended_inputs';O=H/'host/new_catalogue';O.mkdir(exist_ok=True)
sparse_path=R/'41_desert_lineage_function_integration_20260930/host/annotation_audit/desert_full_KO_sparse.tsv';meta_path=R/'41_desert_lineage_function_integration_20260930/host/primary_quality_taxonomy_ANI.tsv';registry_path=R/'40_desert_quality_classification_20260930/host/desert_MAG911_quality_source_registry.tsv'
s=pd.read_csv(sparse_path,sep='\t');allmeta=pd.read_csv(meta_path,sep='\t');m=allmeta[allmeta.dataset_id=='global_desert_PRJNA832417'].copy();reg=pd.read_csv(registry_path,sep='\t');m=m.merge(reg[['assembly','parent_sample_ids','biosample']],left_on='accession',right_on='assembly',validate='one_to_one');assert len(m)==826 and set(m.genome_id)==set(s.genome_id)
old=pd.read_csv(A/'core65_MAG_module_breadth.tsv.gz',sep='\t',usecols=['genome_id']).genome_id.unique();oldrows=allmeta[allmeta.genome_id.isin(old)];assert len(oldrows)==3904
orig_clusters=set(oldrows.ani_cluster_95);anyold_clusters=set(allmeta[allmeta.dataset_id!='global_desert_PRJNA832417'].ani_cluster_95)
m['shares_ANI95_with_original3904']=m.ani_cluster_95.isin(orig_clusters);m['shares_ANI95_with_any_existing']=m.ani_cluster_95.isin(anyold_clusters)
m.to_csv(O/'desert826_endpoint_metadata.tsv',sep='\t',index=False)
summary={'MAGs':len(m),'bacteria':int((m.domain=='Bacteria').sum()),'archaea':int((m.domain=='Archaea').sum()),'sites':m.site_labels.nunique(),'strict_MAGs':int(m.passes_strict_90_5.sum()),'background_flagged_MAGs':int(m.background_flagged.sum()),'MAGs_in_ANI95_shared_original3904':int(m.shares_ANI95_with_original3904.sum()),'MAGs_in_ANI95_shared_any_existing':int(m.shares_ANI95_with_any_existing.sum()),'new_studies':1,'not_external_replication':True,'no_ecological_model_fitted':True}
for branch in ['core65','restored74']:
 p=pd.read_csv(A/f'{branch}_panel.tsv',sep='\t');hits=s[s.KO.isin(p.KO)];assert not hits.duplicated(['genome_id','KO']).any();mat=hits.pivot(index='genome_id',columns='KO',values='protein_count').reindex(index=m.genome_id,columns=p.KO).fillna(0).gt(0).astype(int);assert mat.shape==(826,len(p))
 mat.to_csv(O/f'{branch}_presence.tsv.gz',sep='\t');rows=[]
 for (axis,module),pp in p.groupby(['strategy_axis','module']):
  n=mat[list(pp.KO)].sum(axis=1)
  for gid,v in n.items():rows.append({'genome_id':gid,'strategy_axis':axis,'module':module,'n_panel_KOs':len(pp),'n_detected_KOs':int(v),'module_breadth':v/len(pp)})
 out=pd.DataFrame(rows);assert len(out)==826*9;out.to_csv(O/f'{branch}_module_breadth.tsv.gz',sep='\t',index=False);summary[branch]={'KOs':len(p),'modules':out.module.nunique(),'zero_hit_KOs':mat.columns[mat.sum(axis=0)==0].tolist()}
summary['input_sha256']={str(p.resolve()):hashlib.sha256(p.read_bytes()).hexdigest() for p in [sparse_path,meta_path,registry_path,A/'core65_panel.tsv',A/'restored74_panel.tsv']}
(O/'alignment_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps({k:v for k,v in summary.items() if k!='input_sha256'},indent=2))
