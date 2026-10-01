"""Redraw completed archived phylogeny and breadth data, without refitting."""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.collections import LineCollection

MODULES=['cold_protein_quality','dna_repair_radiation','dormancy_resuscitation','osmotic_desiccation_salt','oxidative_redox','light_energy','sulfur_chemolithotrophy','trace_gas_energy','biofilm_eps_surface']
LABELS=['Cold / protein\nquality','DNA repair /\nprotection','Dormancy /\nstringent','Solute / ion\nhomeostasis','Oxidative /\nredox','Light\nmarkers','Sulfur\nmetabolism','H2 / CO /\ncarbon','Surface /\nEPS']

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo',type=Path,required=True);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 s=a.repo/'data/upstream/panel_repair/figure2_source';d=a.repo/'data/derived'
 means=pd.read_csv(s/'Figure2B_lineage_module_breadth.csv');order=pd.read_csv(s/'Figure2B_tree_tip_order.csv');meta=pd.read_csv(s/'monophyletic_lineage_summary.csv').set_index('lineage_id')
 assert not means.duplicated(['domain','lineage_id','module']).any()
 plt.rcParams.update({'font.family':'Arial','font.size':6.5,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5,'xtick.major.width':.5,'ytick.major.width':.5})
 cmap=LinearSegmentedColormap.from_list('breadth',['#f3f5f5','#c2dddd','#72aeb9','#32738c','#143d62'])
 fig=plt.figure(figsize=(183/25.4,195/25.4),facecolor='white')
 report={};plotted=[];selected=[]
 lefts={'tree':.055,'matrix':.245,'name':.675,'n':.91,'sources':.96}
 widths={'tree':.18,'matrix':.41,'name':.22,'n':.038,'sources':.025}
 for i,(domain,bottom,height) in enumerate([('Bacteria',.315,.525),('Archaea',.163,.072)]):
  ids=order[order.domain.eq(domain)].sort_values('y',ascending=False).lineage_id.tolist();pos={k:j for j,k in enumerate(ids)}
  nodes=pd.read_csv(d/f'{domain.lower()}_tree_nodes.tsv',sep='\t').fillna('');edges=pd.read_csv(d/f'{domain.lower()}_tree_edges.tsv',sep='\t')
  nn=nodes.set_index('node').label.to_dict();depth=nodes.set_index('node').depth.to_dict();children=edges.groupby('parent').child.apply(list).to_dict();ys={}
  assert set(nodes.loc[nodes.label.ne(''),'label'])==set(ids)
  assert len(edges)==len(nodes)-1
  def ypos(n):
   if n not in ys:ys[n]=pos[nn[n]] if nn[n] else np.mean([ypos(c) for c in children[n]])
   return ys[n]
  for n in nodes.node:ypos(n)
  ax={k:fig.add_axes([lefts[k],bottom,widths[k],height]) for k in lefts}
  for v in ax.values():v.set_ylim(len(ids)-.5,-.5)
  segments=[]
  for p,cs in children.items():
   segments.append([(depth[p],min(ys[c] for c in cs)),(depth[p],max(ys[c] for c in cs))])
   for c in cs:segments.append([(depth[p],ys[c]),(depth[c],ys[c])])
  ax['tree'].add_collection(LineCollection(segments,colors='#425761',linewidths=.55))
  xmax=max(depth.values())*1.02
  guides=[[(depth[n],pos[label]),(xmax,pos[label])] for n,label in nn.items() if label]
  ax['tree'].add_collection(LineCollection(guides,colors='#d8dfe2',linewidths=.3,linestyles='dotted'))
  ax['tree'].set(xlim=(-.025*xmax,xmax),yticks=[])
  ax['tree'].spines[['left','right','top','bottom']].set_visible(False);ax['tree'].set_xticks([])
  # The exact branch scale is separate from row alignment guides.
  scale=.5;yy=len(ids)+1.0
  ax['tree'].plot([0,scale],[yy,yy],color='#425761',lw=.8,clip_on=False)
  ax['tree'].text(scale+.06,yy,'0.5',ha='left',va='center',fontsize=5.5,clip_on=False)
  m=means[means.domain.eq(domain)].pivot(index='lineage_id',columns='module',values='mean_module_breadth').loc[ids,MODULES]
  assert m.shape==(len(ids),9) and np.isfinite(m.values).all() and ((m.values>=0)&(m.values<=1)).all()
  im=ax['matrix'].imshow(m.values,aspect='auto',interpolation='nearest',cmap=cmap,vmin=0,vmax=1,extent=(-.5,8.5,len(ids)-.5,-.5))
  ax['matrix'].set(xticks=[],yticks=[]);ax['matrix'].spines[:].set_visible(False)
  for x in [4.5,7.5]:ax['matrix'].axvline(x,color='white',lw=2)
  for k in ['name','n','sources']:ax[k].set_xlim(0,1);ax[k].axis('off')
  chosen=[]
  if i==0:
   for k in meta.loc[ids].sort_values('genomes',ascending=False).index:
    if meta.loc[k,'genomes']>=20 and all(abs(pos[k]-pos[q])>=4 for q in chosen):chosen.append(k)
  else:chosen=ids
  for k in chosen:
   y=pos[k];row=meta.loc[k]
   ax['name'].text(0,y,row.collapsed_name,va='center',fontsize=6.1)
   ax['n'].text(.5,y,str(int(row.genomes)),ha='center',va='center',fontsize=6.1)
   ax['sources'].text(.5,y,str(int(row.source_datasets)),ha='center',va='center',fontsize=6.1,color='#5b686d')
   selected.append({'domain':domain,'lineage_id':k,'name':row.collapsed_name,'row':y,'genomes':int(row.genomes),'sources':int(row.source_datasets)})
  for k in ids:
   for module in MODULES:plotted.append({'domain':domain,'lineage_id':k,'row':pos[k],'module':module,'mean_module_breadth':m.loc[k,module]})
  ytop=bottom+height
  fig.text(.045,ytop+.055 if i==0 else ytop+.031,'a' if i==0 else 'b',weight='bold',fontsize=9)
  fig.text(.08,ytop+.055 if i==0 else ytop+.031,domain,weight='bold',fontsize=9)
  fig.text(.08,ytop+.035 if i==0 else ytop+.013,f'{len(ids)} blocks · {meta.loc[ids,"genomes"].sum():,} MAGs',fontsize=6.5,color='#51636b')
  if i==0:
   for j,label in enumerate(LABELS):
    x=lefts['matrix']+widths['matrix']*(j+.5)/9
    fig.text(x,ytop+.015,label,rotation=55,ha='left',va='bottom',fontsize=6.0)
   for lo,hi,label in [(0,5,'Maintenance'),(5,8,'Energy-associated'),(8,9,'Surface')]:
    x0=lefts['matrix']+widths['matrix']*lo/9;x1=lefts['matrix']+widths['matrix']*hi/9
    fig.add_artist(plt.Line2D([x0+.003,x1-.003],[.948,.948],transform=fig.transFigure,color='#607984',lw=1))
    fig.text((x0+x1)/2,.956,label,ha='center',fontsize=6.4)
  for k,title in [('name','Selected lineages'),('n','MAGs'),('sources','Sources')]:
   fig.text(lefts[k]+(widths[k]/2 if k!='name' else 0),ytop+.011,title,ha='center' if k!='name' else 'left',fontsize=5.8,color='#465d67')
  report[domain]={'blocks':len(ids),'display_genomes':int(meta.loc[ids,'genomes'].sum()),'nodes':len(nodes),'edges':len(edges),'labelled_blocks':len(chosen),'module_values':int(m.size)}
 cax=fig.add_axes([.245,.094,.25,.011]);cb=fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[0,.25,.5,.75,1]);cb.outline.set_visible(False);cb.ax.tick_params(length=2,labelsize=6)
 fig.text(.525,.097,'Mean fraction of specified KOs',fontsize=6.5,va='center')
 fig.text(.055,.016,'Each row is one monophyletic block; terminal branches show archived display representatives.\nColour denotes encoded marker breadth, not activity. Branch scale: substitutions/site.\nSource counts describe catalogue coverage, not independent environmental replicates.',fontsize=6.0,color='#51636b',linespacing=1.55)
 for ext in ['png','pdf','svg']:fig.savefig(a.out/f'figure2_lineage_redesign.{ext}',dpi=600,facecolor='white')
 plt.close(fig)
 frame=pd.DataFrame(plotted);frame.to_csv(a.out/'plotted_matrix.tsv',sep='\t',index=False)
 original=means.set_index(['domain','lineage_id','module']).mean_module_breadth
 candidate=frame.set_index(['domain','lineage_id','module']).mean_module_breadth
 assert np.array_equal(original.sort_index().values,candidate.sort_index().values)
 pd.DataFrame(selected).to_csv(a.out/'display_labels.tsv',sep='\t',index=False)
 report.update({'all_711_values_exact':True,'branch_lengths_and_topology_preserved':True,'cohort':'archived completed core65 primary display blocks','new1129_tree_included':False,'pending607_included':False,'dimensions_mm':[183,195],'reference':'original redesign; WeChat link not accessible','input_hashes':{str(p.relative_to(a.repo)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [s/'Figure2B_lineage_module_breadth.csv',s/'Figure2B_tree_tip_order.csv',s/'monophyletic_lineage_summary.csv',*sorted(d.glob('*tree_nodes.tsv')),*sorted(d.glob('*tree_edges.tsv'))]}})
 (a.out/'validation.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':main()
