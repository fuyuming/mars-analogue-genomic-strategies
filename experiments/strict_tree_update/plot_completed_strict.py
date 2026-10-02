from pathlib import Path
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from Bio import Phylo
H=Path(__file__).resolve().parent;O=H/'host';F=H/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False})
meta=pd.read_csv(O/'strict1243_integrated_metadata.tsv',sep='\t').set_index('tip_id');mp=pd.read_csv(O/'strict_tree_lineage_mapping.tsv',sep='\t');sm=pd.read_csv(O/'strict_lineage_summary.tsv',sep='\t').set_index('lineage_id');breadth=pd.read_csv(O/'strict1243_module_breadth.tsv.gz',sep='\t');breadth=breadth[breadth.branch=='core65'].merge(mp[['tip_id','lineage_id']],on='tip_id',validate='many_to_one')
sources=[('alaska_permafrost_reference','Legacy soil reference'),('australian_basalt_lava_tubes_bay_2025','Australian lava tubes'),('mauna_loa_lava_tube_fishman_2023','Mauna Loa lava tube'),('qaidam_basin_mmag_1773','Qaidam Basin'),('stordalen_mire_2019_hybrid_mags','Stordalen mire'),('atacama_halite_rainfall_prjna484015','Atacama halite'),('atacama_salt_crust_prjna351262','Atacama salt crust'),('atacama_TLT_PRJNA1104199','Atacama TLT'),('atacama_boulder_fields_PRJNA665391','Atacama boulder fields'),('global_desert_PRJNA832417','Global desert catalogue'),('mackay_PRJNA630822','Mackay Glacier soils'),('danakil_accession_catalogue','Danakil catalogue')]
mods=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox','sulfur_chemolithotrophy','trace_gas_energy','light_energy','biofilm_eps_surface'];labels=['Protein quality','DNA repair','Dormancy','Solute / ion','Redox defence','Sulfur energy','Trace-gas marks','Light marks','EPS / surface']
tracks=[]
for marker,domain in [('bac120','Bacteria'),('ar53','Archaea')]:
 tree=Phylo.read(O/f'{marker}_collapsed_display.newick','newick');leaves=tree.get_terminals();order=[t.name for t in leaves];N=len(order);height=7.7 if N>10 else 4.5
 fig=plt.figure(figsize=(7.1,height));bottom=.26 if N>10 else .51;top=.91;span=top-bottom
 ax=fig.add_axes([.025,bottom,.405,span]);an=fig.add_axes([.444,bottom,.034,span]);asc=fig.add_axes([.495,bottom,.211,span]);af=fig.add_axes([.735,bottom,.24,span]);depth=tree.depths();maxx=max(depth[t] for t in leaves);yy={t:i for i,t in enumerate(leaves)}
 seg=[]
 for c in tree.find_clades(order='postorder'):
  if c.is_terminal():continue
  yy[c]=sum(yy[k] for k in c.clades)/len(c.clades)
  seg.append([(depth[c],min(yy[k] for k in c.clades)),(depth[c],max(yy[k] for k in c.clades))])
  for k in c.clades:seg.append([(depth[c],yy[k]),(depth[k],yy[k])])
 ax.add_collection(LineCollection(seg,colors='#67757D',linewidths=.55))
 for t in leaves:
  label=sm.loc[t.name,'collapsed_name'];label=('c. ' if sm.loc[t.name,'collapsed_rank']=='class' else '')+label
  ax.plot([depth[t],maxx*1.02],[yy[t],yy[t]],color='#D6DADD',lw=.35,ls=':');ax.text(maxx*1.04,yy[t],label,va='center',fontsize=6.6 if N>10 else 7)
 ax.set_xlim(-maxx*.02,maxx*2.02);ax.set_ylim(N-.5,-.5);ax.axis('off');scale=.2;ax.plot([0,scale],[N+.2,N+.2],color='#37454D',lw=.9,clip_on=False);ax.text(0,N+.65,'0.2 substitutions/site',ha='left',va='top',fontsize=6,clip_on=False)
 for i,lid in enumerate(order):an.text(.5,i,str(sm.loc[lid,'genomes']),ha='center',va='center',fontsize=6.8)
 an.set(xlim=(0,1),ylim=(N-.5,-.5));an.axis('off');an.set_title('MAGs',fontsize=7,pad=8)
 joined=mp[mp.domain==domain].merge(meta[['dataset_id']],left_on='tip_id',right_index=True,validate='one_to_one');ct=pd.crosstab(joined.lineage_id,joined.dataset_id).reindex(index=order,columns=[s[0] for s in sources],fill_value=0);sc=ct.div(ct.sum(axis=1),axis=0);fm=breadth.groupby(['lineage_id','module']).breadth.mean().unstack().loc[order,mods]
 for i,lid in enumerate(order):
  for j,(src,_) in enumerate(sources):tracks.append(dict(domain=domain,lineage_id=lid,track='source_fraction',feature=src,n=int(ct.loc[lid,src]),denominator=int(ct.loc[lid].sum()),value=float(sc.loc[lid,src])))
  for mod in mods:tracks.append(dict(domain=domain,lineage_id=lid,track='mean_marker_breadth',feature=mod,n=int(sm.loc[lid,'genomes']),denominator=int(breadth.loc[(breadth.lineage_id==lid)&(breadth.module==mod),'denominator'].iloc[0]),value=float(fm.loc[lid,mod])))
 im=asc.imshow(sc.to_numpy(),aspect='auto',vmin=0,vmax=1,cmap='Blues',interpolation='nearest');im2=af.imshow(fm.to_numpy(),aspect='auto',vmin=0,vmax=1,cmap='YlGnBu',interpolation='nearest')
 for a,ls in [(asc,[f'S{i+1}' for i in range(12)]),(af,labels)]:
  a.set_yticks([]);a.set_xticks(range(len(ls)),ls,rotation=90,fontsize=6);a.tick_params(length=0);a.set_ylim(N-.5,-.5)
  for sp in a.spines.values():sp.set_visible(False)
  a.set_xticks(np.arange(-.5,len(ls),1),minor=True);a.set_yticks(np.arange(-.5,N,1),minor=True);a.grid(which='minor',color='white',lw=.3);a.tick_params(which='minor',length=0)
 asc.set_title('Catalogue composition',fontsize=8,pad=8);af.set_title('Encoded marker breadth',fontsize=8,pad=8);af.axvline(4.5,color='#223F4A',lw=.85);af.axvline(7.5,color='#223F4A',lw=.6)
 fig.text(.025,.965,domain+' | expanded strict cohort',weight='bold',fontsize=11);fig.text(.025,.933,f'{int(sm.loc[order,"genomes"].sum()):,} MAGs · {N} sequence-tree groups · core65 panel',fontsize=8,color='#50616A')
 legend_y=.108 if N>10 else .20
 for i,(_,name) in enumerate(sources):
  col=i//4;row=i%4;fig.text(.025+col*.325,legend_y-row*(.018 if N>10 else .038),f'S{i+1}  {name}',fontsize=6.4)
 cy=.153 if N>10 else .30
 for xpos,obj,title in [(.49,im,'Within-group fraction'),(.735,im2,'Mean fraction detected')]:
  cax=fig.add_axes([xpos,cy,.20,.009 if N>10 else .018]);cb=fig.colorbar(obj,cax=cax,orientation='horizontal',ticks=[0,.5,1]);cb.ax.tick_params(labelsize=6,length=2);cb.set_label(title,fontsize=6,labelpad=1);cb.outline.set_visible(False)
 fig.text(.025,.015,'Sequence phylogram scaffold; rooting is for display. Catalogue ≠ independent habitat replicate. Breadth ≠ activity.',fontsize=6,color='#50616A')
 for ext in ['svg','pdf','png']:fig.savefig(F/f'strict_{marker}_tree_function.{ext}',dpi=250)
 plt.close(fig)
pd.DataFrame(tracks).to_csv(O/'tree_figure_source_data.tsv',sep='\t',index=False)
# Paired uncertainty: dot is bootstrap median, not an ecological posterior mean.
s=pd.read_csv(O/'strict_source_panel_sensitivity.tsv',sep='\t');s['label']=s.apply(lambda r:f"{r.branch} · "+('all maintenance' if r.panel=='full' else 'without solute/ion')+' · '+('lineage weights' if r.scheme=='lineage_BB' else 'lineage + catalogue weights'),axis=1);fig,ax=plt.subplots(figsize=(7.1,3.3));fig.subplots_adjust(left=.51,right=.96,bottom=.20,top=.84)
for i,r in s.iterrows():
 color='#226D82' if r.branch=='core65' else '#8A6337';ax.plot([100*r.q025,100*r.q975],[i,i],color=color,lw=1.6);ax.scatter(100*r.delta_median,i,color=color,s=22,zorder=3)
ax.axvline(0,color='#7D858B',lw=.8,ls='--');ax.set_yticks(range(len(s)),s.label,fontsize=6.7);ax.invert_yaxis();ax.set_xlabel('Maintenance − energy source partial R² (percentage points)',fontsize=8);ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='y',length=0);ax.grid(axis='x',color='#E7EAEC',lw=.5);fig.text(.045,.93,'Strict-cohort source contrast',fontsize=11,weight='bold');fig.text(.045,.87,'1,154 bacteria · 23 shared tree groups · 12 catalogues',fontsize=8);fig.text(.045,.03,'Dots: paired Bayesian-bootstrap medians; bars: 95% weight-distribution intervals. Conditional on observed catalogues.',fontsize=6.3)
for ext in ['svg','pdf','png']:fig.savefig(F/f'strict_paired_source_contrast.{ext}',dpi=250)
plt.close(fig)
