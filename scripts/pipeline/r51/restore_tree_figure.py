from pathlib import Path
import pandas as pd,numpy as np,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;O=H/'host';F=H/'figures';D=H/'manuscript/assets';A=H.parents[1].parent/'code/result_raw/KO_panel_repair_20260926/word_revision/delivery'
for ext in ['png','pdf','svg']:shutil.copy2(F/f'figure2_lineage_structure.{ext}',F/f'figureS3_null_details.{ext}')
shutil.copy2(F/'figureS3_null_details.png',D/'figureS3_null_details.png')
shutil.copy2(A/'Figures_corrected/Figure5_core65_corrected.png',D/'figureS4_independent_occurrence.png')
plt.rcParams.update({'font.family':'Arial','font.size':12,'axes.titlesize':13,'axes.labelsize':12,'xtick.labelsize':11,'ytick.labelsize':11,'pdf.fonttype':42,'svg.fonttype':'none'})
means=pd.read_csv(A/'Source_Data/core65/Figure2B_lineage_module_breadth.csv');order=pd.read_csv(A/'Source_Data/core65/Figure2B_tree_tip_order.csv');meta=pd.read_csv(A/'Source_Data/core65/monophyletic_lineage_summary.csv').set_index('lineage_id')
modules=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox','light_energy','sulfur_chemolithotrophy','trace_gas_energy','biofilm_eps_surface']
labels=['Cold / protein\nquality','DNA repair /\nprotection','Dormancy /\nstringent','Solute / ion\nhomeostasis','Oxidative /\nredox','Light-\nassociated','Sulfur\nmetabolism','H2 / CO /\ncarbon','Surface /\nEPS']
fig=plt.figure(figsize=(10,7.6),layout='constrained');gs=fig.add_gridspec(2,3,height_ratios=[4.6,1.2],width_ratios=[1.3,1.6,5.2],hspace=.14,wspace=.025)
matrix_out=[];display=[]
for i,dom in enumerate(['Bacteria','Archaea']):
 ax=fig.add_subplot(gs[i,0]);la=fig.add_subplot(gs[i,1]);hm=fig.add_subplot(gs[i,2]);ids=order[order.domain==dom].sort_values('y',ascending=False).lineage_id.tolist();pos={k:j for j,k in enumerate(ids)}
 nodes=pd.read_csv(O/f'{dom.lower()}_tree_nodes.tsv',sep='\t').fillna('');edges=pd.read_csv(O/f'{dom.lower()}_tree_edges.tsv',sep='\t');depth=nodes.set_index('node').depth.to_dict();nn=nodes.set_index('node').label.to_dict();child=edges.groupby('parent').child.apply(list).to_dict();ys={}
 def calc(n):
  if n in ys:return ys[n]
  ys[n]=pos[nn[n]] if nn[n] else float(np.mean([calc(c) for c in child[n]]));return ys[n]
 for n in nodes.node:calc(n)
 for p,cs in child.items():
  ax.plot([depth[p]]*2,[min(ys[c] for c in cs),max(ys[c] for c in cs)],color='#69777E',lw=.6)
  for c in cs:ax.plot([depth[p],depth[c]],[ys[c]]*2,color='#69777E',lw=.6)
 ax.set(ylim=(len(ids)-.5,-.5),xlabel='Branch length',yticks=[]);ax.spines[['top','right','left']].set_visible(False);ax.tick_params(axis='x',labelsize=9)
 ax.set_title(('a' if i==0 else 'b')+'  '+dom+'\n'+str(len(ids))+' lineage blocks',loc='left',fontsize=13,fontweight='bold')
 la.set(xlim=(0,1),ylim=(len(ids)-.5,-.5));la.axis('off');chosen=[]
 if i==0:
  for k in meta.loc[ids].sort_values('genomes',ascending=False).index:
   if meta.loc[k,'genomes']>=20 and all(abs(pos[k]-pos[x])>=4 for x in chosen):chosen.append(k)
 else:chosen=ids
 for k in chosen:la.text(.98,pos[k],str(meta.loc[k,'collapsed_name'])+' ('+str(meta.loc[k,'genomes'])+')',ha='right',va='center',fontsize=10.5);display.append(dict(domain=dom,lineage_id=k,name=meta.loc[k,'collapsed_name'],genomes=meta.loc[k,'genomes'],row=pos[k]))
 m=means[means.domain==dom].pivot(index='lineage_id',columns='module',values='mean_module_breadth').loc[ids,modules];matrix_out.append(m.assign(domain=dom));im=hm.imshow(m,aspect='auto',vmin=0,vmax=1,cmap='YlGnBu',interpolation='nearest');hm.set(yticks=[],xticks=range(9));hm.set_xticklabels(labels if i==1 else [],rotation=45,ha='right');hm.tick_params(axis='both',length=0);hm.spines[:].set_visible(False)
 if i==0:
  hm.set_title('Mean encoded marker breadth by lineage',loc='left',fontsize=13,fontweight='bold')
  for x in [4.5,7.5]:hm.axvline(x,color='white',lw=2)
 else:
  for x in [4.5,7.5]:hm.axvline(x,color='white',lw=2)
fig.colorbar(im,ax=fig.axes,location='right',shrink=.65,pad=.02,label='Mean fraction of specified KOs')
for ext in ['png','pdf','svg']:fig.savefig(F/f'figure2_lineage_mosaic.{ext}',dpi=240,bbox_inches='tight',facecolor='white')
shutil.copy2(F/'figure2_lineage_mosaic.png',D/'figure2_lineage_mosaic.png');pd.concat(matrix_out).to_csv(O/'figure2_lineage_mosaic.tsv',sep='\t');pd.DataFrame(display).to_csv(O/'figure2_display_labels.tsv',sep='\t',index=False)
print('Restored original archived tree/module evidence; no tree reconstruction or statistical rerun.')
