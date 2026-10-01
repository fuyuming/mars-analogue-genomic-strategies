from pathlib import Path
import pandas as pd,numpy as np,json,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow,Patch
H=Path(__file__).resolve().parent;R=H.parent;A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926';O=H/'host';F=H/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Arial','font.size':11,'axes.titlesize':12,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
teal='#176B80';orange='#BB6C32';grey='#B3B9BC'
sources=['qaidam_basin_mmag_1773','alaska_permafrost_reference','stordalen_mire_2019_hybrid_mags','atacama_salt_crust_prjna351262','atacama_halite_rainfall_prjna484015','mauna_loa_lava_tube_fishman_2023','australian_basalt_lava_tubes_bay_2025']
labels=['Qaidam','Alaska','Stordalen','Atacama crust','Atacama halite','Mauna Loa','Australian caves']
sd=dict(zip(sources,labels));systems=['KdpABC','KdpDE','ProVWX','EctABC','EctD','BetAB']
module_names={'cold_protein_quality':'Cold / protein quality','dna_repair_radiation':'DNA repair / protection','dormancy_resuscitation':'Dormancy / stringent','osmotic_desiccation_salt':'Solute / ion homeostasis','oxidative_redox':'Oxidative / redox','light_energy':'Light-associated','sulfur_chemolithotrophy':'Sulfur metabolism','trace_gas_energy':'H2 / CO / carbon markers','biofilm_eps_surface':'Surface / EPS'}
def title(ax,letter,s):ax.set_title(letter+'  '+s,loc='left',fontweight='bold',pad=12)
def save(fig,name):
 for ext in ['pdf','svg','png']:fig.savefig(F/f'{name}.{ext}',dpi=230,bbox_inches='tight',facecolor='white')
 plt.close(fig);shutil.copy2(F/f'{name}.png',H/'manuscript/assets'/f'{name}.png')
def table(d,n):d.to_csv(O/n,sep='\t',index=False)

b=pd.read_csv(A/'amended_inputs/core65_MAG_module_breadth.tsv.gz',sep='\t');meta=b.drop_duplicates('genome_id');rep=meta[meta.species_representative].copy()
counts=meta.groupby('dataset_id').agg(primary=('genome_id','size'),strict=('passes_strict_90_5','sum')).reindex(sources);counts['ANI95']=rep.groupby('dataset_id').size();table(counts.reset_index(),'figure1_cohort_counts.tsv')
top=list(rep.phylum.value_counts().head(5).index)+['Halobacteriota'];rep['display_phylum']=np.where(rep.phylum.isin(top),rep.phylum,'Other / unclassified');ct=pd.crosstab(rep.dataset_id,rep.display_phylum).reindex(sources);ct=ct.reindex(columns=list(top)+['Other / unclassified']);frac=ct.div(ct.sum(1),axis=0)
table(ct.reset_index(),'figure1_taxonomic_counts.tsv')
co=pd.read_csv(R/'48_MAG_portfolio_inference_20261001/host/core65_primary_cohort.tsv',sep='\t');means=b[b.genome_id.isin(co.genome_id)].groupby(['dataset_id','module']).module_breadth.mean().unstack().reindex(index=sources,columns=list(module_names));table(means.reset_index(),'figure1_module_means.tsv')
fig=plt.figure(figsize=(10,7.5),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.1]);a=fig.add_subplot(gs[0,0]);c=fig.add_subplot(gs[0,1]);e=fig.add_subplot(gs[1,:]);y=np.arange(7)
a.barh(y,counts.primary,color='#DCE3E5',label='Primary MAGs');a.barh(y,counts.strict,color=teal,label='Strict subset');a.plot(counts.ANI95,y,'D',ms=4.5,color=orange,label='ANI95 representatives')
for i,v in enumerate(counts.primary):a.text(v+20,i,str(v),va='center',fontsize=9)
a.set(yticks=y,yticklabels=labels,xlabel='Genomes',xlim=(0,1900));a.invert_yaxis();a.legend(frameon=False,fontsize=9,loc='center right',bbox_to_anchor=(1,.40));title(a,'a','Catalogue representation')
left=np.zeros(7);colors=['#176B80','#5F9F9D','#DEA158','#A67EAB','#7C8C53','#887455','#D9D9D9']
for col,colr in zip(frac,colors):c.barh(y,frac[col],left=left,color=colr,label=col);left+=frac[col].values
c.set(yticks=y,yticklabels=[],xlabel='Fraction of ANI95 representatives',xlim=(0,1));c.invert_yaxis();title(c,'b','Taxonomic composition');c.legend(frameon=False,fontsize=9,bbox_to_anchor=(1.01,.98),loc='upper left')
im=e.imshow(means.values.T,aspect='auto',vmin=0,vmax=1,cmap='YlGnBu');e.set(xticks=y,xticklabels=labels,yticks=range(9),yticklabels=[module_names[m] for m in means]);e.tick_params(axis='x',labelrotation=18)
for i in range(9):
 for j in range(7):e.text(j,i,('<0.01' if 0<means.iloc[j,i]<.005 else f'{means.iloc[j,i]:.2f}'),ha='center',va='center',fontsize=9,color='white' if means.iloc[j,i]>.6 else 'black')
title(e,'c','Encoded module breadth in the bacterial comparison');fig.colorbar(im,ax=e,shrink=.8,label='Mean fraction of specified KOs');save(fig,'figure1_cohort_strategies')

r=pd.read_csv(A/'recurrence_core65/KO_recurrence_null_results.csv');r=r[(r.stratum=='primary_ANI95_representatives')&r.variable_trait];rn=r[r.null_model=='within_source_prevalence_preserving'];table(rn,'figure2_KO_null_data.tsv')
fig=plt.figure(figsize=(10,5.6),layout='constrained');gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.12])
for j,dom in enumerate(['Bacteria','Archaea']):
 ax=fig.add_subplot(gs[0,j]);d=rn[rn.domain==dom];lim=0
 for axis,col,mark,lbl in [('cell_maintenance',teal,'o','Maintenance'),('energy_acquisition',orange,'^','Energy'),('surface_retention',grey,'s','Surface')]:
  g=d[d.strategy_axis==axis];x=g.null_mean/g.n_tips;y=g.observed_score/g.n_tips;ax.scatter(x,y,s=32,color=col,marker=mark,alpha=.75,label=lbl);lim=max(lim,x.max(),y.max())
 ax.plot([0,lim*1.06],[0,lim*1.06],'--',color='.6',lw=1);ax.set(xlabel='Null mean transitions / tip',ylabel='Observed transitions / tip',xlim=(-.006,lim*1.06),ylim=(-.006,lim*1.06));title(ax,chr(97+j),dom)
 if j==0:ax.legend(frameon=False,fontsize=9,loc='upper left')
ax=fig.add_subplot(gs[0,2]);summary=[]
for dom in ['Bacteria','Archaea']:
 for axis in ['cell_maintenance','energy_acquisition']:
  for null in ['global_prevalence_preserving','within_source_prevalence_preserving']:
   d=r[(r.domain==dom)&(r.strategy_axis==axis)&(r.null_model==null)];summary.append(dict(domain=dom,axis=axis,null=null,n=len(d),clustered=int((d.p_less_bh<.05).sum())))
ss=pd.DataFrame(summary);table(ss,'figure2_clustered_counts.tsv')
for i,(dom,axis) in enumerate([(d,a) for d in ['Bacteria','Archaea'] for a in ['cell_maintenance','energy_acquisition']]):
 for j,null in enumerate(['global_prevalence_preserving','within_source_prevalence_preserving']):
  row=ss[(ss.domain==dom)&(ss.axis==axis)&(ss['null']==null)].iloc[0];y=i+(j-.5)*.3;ax.barh(y,row.clustered/row.n,height=.26,color=teal if j else '#BBCDD2');ax.text(.02,y,f'{row.clustered}/{row.n}',va='center',fontsize=9,color='white' if j else 'black')
ax.set(yticks=range(4),yticklabels=['Bacteria\nMaintenance','Bacteria\nEnergy','Archaea\nMaintenance','Archaea\nEnergy'],xlabel='Fraction clustered (BH q < 0.05)',xlim=(0,1.02));ax.invert_yaxis();title(ax,'c','Clustering persists within source');ax.legend(handles=[Patch(color='#BBCDD2',label='Global null'),Patch(color=teal,label='Within-source null')],frameon=False,fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.16));save(fig,'figure2_lineage_structure')

# Preserve the exact established paired-inference data figure, including caveats.
for ext in ['png','pdf']:shutil.copy2(R/f'48_MAG_portfolio_inference_20261001/host/portfolio_uncertainty_diagnostic.{ext}',F/f'figure3_portfolio_uncertainty.{ext}')
shutil.copy2(F/'figure3_portfolio_uncertainty.png',H/'manuscript/assets/figure3_portfolio_uncertainty.png')

old=pd.read_csv(R/'49_system_configuration_20261001/host/original_primary_system_states.tsv',sep='\t');m=old.groupby('dataset_id')[systems].agg(lambda x:(x=='panel_complete').mean()).reindex(sources);table(m.reset_index(),'figure4_source_codetection.tsv')
fig=plt.figure(figsize=(10,6.3),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.25,1],height_ratios=[1.2,1]);a=fig.add_subplot(gs[0,:]);c=fig.add_subplot(gs[1,0]);e=fig.add_subplot(gs[1,1])
im=a.imshow(m.values.T*100,aspect='auto',vmin=0,vmax=100,cmap='YlGnBu');a.set(xticks=range(7),xticklabels=[label+'\n(n='+str(int((old.dataset_id==src).sum()))+')' for src,label in zip(sources,labels)],yticks=range(6),yticklabels=systems)
for i in range(6):
 for j in range(7):a.text(j,i,f'{m.iloc[j,i]*100:.1f}',ha='center',va='center',color='white' if m.iloc[j,i]>.6 else 'black',fontsize=10)
title(a,'a','Specific marker combinations vary across original sources');fig.colorbar(im,ax=a,shrink=.9,label='Genomes with all specified markers (%)')
assoc=pd.read_csv(R/'49_system_configuration_20261001/host/lineage_group_associations.tsv',sep='\t');aa=assoc[(assoc.cohort=='original_primary')&(assoc['rank']=='family')].set_index('endpoint').loc[[s+'_codetected' for s in systems]];table(aa.reset_index(),'figure4_family_associations.tsv')
c.plot(aa.source_or_site_partial_r2*100,range(6),'o',color=teal,label='Raw partial R²');c.plot(aa.df_adjusted_partial_r2*100,range(6),'x',color=orange,label='Degrees-of-freedom adjusted');c.axvline(0,color='.75',lw=.7);c.set(yticks=range(6),yticklabels=systems,xlabel='Source partial R² (%)',xlim=(-2,25));c.invert_yaxis();title(c,'b','Within 30 supported families');c.legend(frameon=False,fontsize=9,loc='lower right')
cx=pd.read_csv(O/'system_context_summary.tsv',sep='\t');cc=cx[(cx.cohort=='original_strict')&(cx.metric=='order_compact')].set_index('system').loc[[s for s in systems if s!='EctD']]
e.barh(range(5),cc.fraction,color=teal,height=.55)
for i,row in enumerate(cc.itertuples()):e.text(row.fraction+.025,i,f'{row.n}/{row.denominator}',va='center',fontsize=10)
e.set(yticks=range(5),yticklabels=cc.index,xlabel='Fraction among co-detected genomes',xlim=(0,1.18));e.invert_yaxis();title(e,'c','Original CDS-order proximity proxy');save(fig,'figure4_molecular_configurations')

fig=plt.figure(figsize=(10.5,7),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.1],width_ratios=[1,1.18]);a=fig.add_subplot(gs[0,0]);c=fig.add_subplot(gs[0,1]);e=fig.add_subplot(gs[1,0]);g=fig.add_subplot(gs[1,1])
ss=pd.read_csv(R/'49_system_configuration_20261001/host/system_state_summary.tsv',sep='\t');ss=ss[ss.cohort=='new_desert_strict'].set_index('system').loc[[s for s in systems if s!='EctD']];x=np.arange(5)
a.barh(x,ss.panel_complete,color=teal,label='All specified markers');a.barh(x,ss.partial,left=ss.panel_complete,color='#D4D8D9',label='Partial marker set')
for i,row in enumerate(ss.itertuples()):a.text(row.panel_complete+row.partial+1,i,f'{row.panel_complete}/{row.panel_complete+row.partial}',va='center',fontsize=9)
a.set(yticks=x,yticklabels=ss.index,xlabel='Genomes with any marker; labels = all / any',xlim=(0,95));a.invert_yaxis();title(a,'a','Strict dryland marker sets (n = 154)');a.legend(frameon=False,fontsize=9,loc='center right',bbox_to_anchor=(1,.39))
cc=cx[(cx.cohort=='new_desert_strict')&(cx.metric=='compact_10kb')].set_index('system').loc[ss.index]
for i,row in enumerate(cc.itertuples()):c.plot([row.conditional_q025,row.conditional_q975],[i,i],color=teal,lw=2);c.plot(row.fraction,i,'o',color=teal);c.text(1.035,i,f'{row.n}/{row.denominator}',va='center',fontsize=10)
c.set(yticks=x,yticklabels=ss.index,xlabel='Fraction among co-detected genomes',xlim=(0,1.26));c.invert_yaxis();title(c,'b','Same strand and contig within 10 kb')
support=pd.read_csv(R/'49_system_configuration_20261001/host/lineage_group_support.tsv',sep='\t');ranks=['class','order','family'];y=np.arange(3)
for j,cohort in enumerate(['new_desert_primary','new_desert_strict']):
 s=support[support.cohort==cohort].set_index('rank').loc[ranks];e.barh(y+(j-.5)*.3,s.n,height=.27,color=teal if j==0 else orange,label='Primary' if j==0 else 'Strict')
 for i,n in enumerate(s.n):e.text(n+8,i+(j-.5)*.3,str(n),va='center',fontsize=10)
e.set(yticks=y,yticklabels=['Class','Order','Family'],xlabel='Representatives in supported lineages',xlim=(0,675));e.invert_yaxis();title(e,'c','Cross-site lineage support');e.legend(frameon=False,fontsize=9,loc='lower right')
loci=pd.read_csv(O/'strict_compact_locus_candidates.tsv',sep='\t');selected=[];cols=['#176B80','#DEA158','#A67EAB']
genes={'K01546':'kdpA','K01547':'kdpB','K01548':'kdpC','K02000':'proV','K02001':'proW','K02002':'proX','K06718':'ectA','K00836':'ectB','K06720':'ectC','K00108':'betA','K00130':'betB'}
for i,system in enumerate(['KdpABC','ProVWX','EctABC','BetAB']):
 candidates=loci[loci.system==system];gid=sorted(candidates.genome_id.unique())[0];d=candidates[candidates.genome_id==gid].sort_values('start');selected.append(d);offset=d.start.min()
 for j,row in enumerate(d.itertuples()):
  st=(row.start-offset)/1000;en=(row.end-offset+1)/1000;sign=1 if row.strand=='+' else -1;beg=st if sign==1 else en
  g.add_patch(FancyArrow(beg,i,sign*(en-st),0,width=.17,head_width=.26,head_length=min(.15,(en-st)*.35),length_includes_head=True,facecolor=cols[j%3],edgecolor='none'))
  g.text((st+en)/2,i-.22,genes[row.KO],ha='center',fontsize=9)
 g.text(0,i+.31,gid.replace('desert2025__','')+' · '+str(d.source_or_site.iloc[0]),fontsize=8.5,color='.35')
g.set(yticks=range(4),yticklabels=['KdpABC','ProVWX','EctABC','BetAB'],xlabel='Distance from first marked CDS (kb)',xlim=(-.05,6.4),ylim=(3.7,-.6));title(g,'d','Illustrative observed loci');g.spines['left'].set_visible(False);g.tick_params(axis='y',length=0)
table(pd.concat(selected),'figure5_displayed_loci.tsv');save(fig,'figure5_dryland_context')

for src,dst in [('figure1_framework','figureS1_framework'),('figureS1_bayesian_profiles','figureS2_bayesian_profiles')]:shutil.copy2(R/f'50_analogue_scope_manuscript_20261001/manuscript/assets/{src}.png',H/f'manuscript/assets/{dst}.png')
print('Five data figures generated; conceptual framework and existing Bayesian profile retained as supplements.')
