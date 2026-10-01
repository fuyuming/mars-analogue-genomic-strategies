from pathlib import Path
import pandas as pd,numpy as np,re,json,itertools,hashlib
H=Path(__file__).resolve().parent;B=H.parent;O=H/'host';O.mkdir(exist_ok=True)
p=B/'37_desert_inclusion_20260930/derived/desert_extension_registry.tsv';t=pd.read_csv(p,sep='\t');t=t[t.sample_role.eq('field_desert_soil')].copy()
p49=B/'49_system_configuration_20261001/host'
a=pd.read_csv(p49/'new_desert_primary_independent_representatives.tsv',sep='\t');s=pd.read_csv(p49/'new_desert_strict_independent_representatives.tsv',sep='\t')
e=pd.read_csv(B/'39_iterative_mechanism_analysis_20260930/host/desert_assembly_parent_edges.tsv',sep='\t')
focus={'Atacama Desert','McMurdo Dry Valleys','Mount Seuss','Mojave Desert','Namib Desert','Negev Desert'}
def ll(x):
 m=re.findall(r'([\d.]+)\s*([NSEW])',str(x));assert len(m)==2,x
 return tuple(float(v)*(-1 if k in 'SW' else 1) for v,k in m)
def km(a,b):
 x,y,u,v=np.radians([*a,*b]);q=np.sin((u-x)/2)**2+np.cos(x)*np.cos(u)*np.sin((v-y)/2)**2
 return float(6371.0088*2*np.arcsin(np.sqrt(np.clip(q,0,1))))
rows=[]
for n,g in t.groupby('site_label_from_title'):
 pts=[ll(x) for x in g.lat_lon];d=[km(x,y) for x,y in itertools.combinations(pts,2)]
 rows.append(dict(site_label=n,regional_literature_focus=n in focus,scope_status='regional_candidate_pending_literature_microhabitat_audit' if n in focus else 'broader_dryland_context_unresolved_analogue_scope',field_samples=len(g),biosamples=';'.join(g.sample_accession),deposited_locations=';'.join(g.location.unique()),deposited_coordinates=';'.join(g.lat_lon.unique()),unique_deposited_coordinates=g.lat_lon.nunique(),max_parent_separation_km=max(d),depths=';'.join(g.depth.unique()),dates=';'.join(g.collection_date.astype(str).unique()),deposited_biomes=';'.join(g.biome_as_deposited.unique()),primary_representatives=int(a.site_labels.eq(n).sum()),strict_representatives=int(s.site_labels.eq(n).sum()),interpretation='deposited metadata geometry; no coordinate correction or ecological independence assumed'))
r=pd.DataFrame(rows);r.to_csv(O/'site_scope_sampling.tsv',sep='\t',index=False)
# Validate actual parent edges for each selected assembly, not inferred label membership.
edges=[]
for cohort,x in [('primary',a),('strict',s)]:
 for _,m in x.iterrows():
  parents=e[e.assembly.eq(m.assembly)].parent_biosample.drop_duplicates();g=t[t.sample_accession.isin(parents)]
  assert len(g)==len(parents) and len(g)>0,(m.assembly,len(g),len(parents))
  assert set(g.site_label_from_title)=={m.site_labels},m.assembly
  coords=[ll(v) for v in g.lat_lon];diam=max([km(u,v) for u,v in itertools.combinations(coords,2)] or [0])
  edges.append(dict(cohort=cohort,genome_id=m.genome_id,assembly=m.assembly,site=m.site_labels,parent_count=len(g),parent_diameter_km=diam,parent_dates=';'.join(g.collection_date.astype(str).unique()),parent_biosamples=';'.join(g.sample_accession),regional_focus=m.site_labels in focus))
v=pd.DataFrame(edges);v.to_csv(O/'selected_assembly_sampling_geometry.tsv',sep='\t',index=False)
# Gate depends on metadata only; no KO or system table read.
out=[]
for cohort,x in [('primary',a),('strict',s)]:
 for scope,z in [('all_25',x),('six_regional_candidates',x[x.site_labels.isin(focus)]),('other_19_context',x[~x.site_labels.isin(focus)])]:
  for rank in ['class','order','family','genus']:
   zz=z[z[rank].notna() & ~z[rank].astype(str).str.fullmatch(r'(?:[a-z]__)?(?:unknown|unclassified|nan|NA)?',case=False)].copy()
   tab=pd.crosstab(zz[rank],zz.site_labels);eligible=tab.index[(tab.ge(5)).sum(axis=1).ge(2)]
   q=zz[zz[rank].isin(eligible)]
   out.append(dict(cohort=cohort,scope=scope,rank=rank,total_representatives=len(z),total_sites=z.site_labels.nunique(),support_representatives=len(q),support_lineages=len(eligible),support_sites=q.site_labels.nunique()))
u=pd.DataFrame(out);u.to_csv(O/'scope_lineage_support.tsv',sep='\t',index=False)
summary={'field_samples':len(t),'labels':t.site_label_from_title.nunique(),'focus_labels':len(focus),'focus_primary':int(a.site_labels.isin(focus).sum()),'focus_strict':int(s.site_labels.isin(focus).sum()),'all_selected_assemblies_parents_verified':True,'assemblies_with_three_parents':int(v[v.cohort.eq('primary')].parent_count.eq(3).sum()),'labels_over_5km':r.loc[r.max_parent_separation_km.gt(5),'site_label'].tolist(),'labels_over_1km':r.loc[r.max_parent_separation_km.gt(1),'site_label'].tolist(),'no_functional_response_loaded':True}
(O/'sampling_audit_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));print(u[u.scope.eq('six_regional_candidates')].to_string(index=False));print(r[['site_label','max_parent_separation_km']].sort_values('max_parent_separation_km',ascending=False).head(8).to_string(index=False))
