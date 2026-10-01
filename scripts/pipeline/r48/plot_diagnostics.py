from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent/'host';d=pd.read_csv(H/'paired_uncertainty.tsv',sep='\t');m=pd.read_csv(H/'module_decomposition.tsv',sep='\t');eq=pd.read_csv(H/'equal_submodule_uncertainty.tsv',sep='\t');points=pd.read_csv(H/'point_reproduction.tsv',sep='\t');influence=pd.read_csv(H/'influence_diagnostics.tsv',sep='\t')
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(13.7,4.6),gridspec_kw={'width_ratios':[1.05,1.1,1.15]})
labels=[]
for i,(branch,strict) in enumerate([('core65',False),('core65',True),('restored74',False),('restored74',True)]):
 labels.append(branch+(' strict' if strict else ' primary'))
 for scheme,off,color in [('lineage_BB',.12,'#1b6a82'),('lineage_source_product_sensitivity',-.12,'#b26b36')]:
  z=d[(d.branch==branch)&(d.strict==strict)&(d.scheme==scheme)&(d.comparison=='full')].iloc[0];ax[0].plot([z.delta_q025*100,z.delta_q975*100],[i+off]*2,color=color,lw=2);ax[0].plot(z.delta_median*100,i+off,'o',color=color,ms=4,label=('Lineage weights' if off>0 else 'Lineage × source weights') if i==0 else None)
ax[0].set(yticks=range(4),yticklabels=labels,xlabel='Maintenance − energy partial R² (pp)',title='A  Paired uncertainty');ax[0].set_ylim(4.3,-.5);ax[0].legend(frameon=False,fontsize=8,loc='lower right');ax[0].axvline(0,c='.6',lw=.8,ls='--')
z=m[(m.branch=='core65') & ~m.strict];names=['Cold-associated','DNA repair/protection','Sporulation/stringent','Solute/ion homeostasis','Oxidative/redox','Sulfur metabolism','H₂/CO/carbon markers'];ax[1].barh(range(7),100*z.source_partial_r2,color=['#1b6a82']*5+['#b26b36']*2);ax[1].set(yticks=range(7),yticklabels=names,xlabel='Module source partial R² (%)',title='B  Module decomposition');ax[1].invert_yaxis()
for i,(comp,label) in enumerate([('full','Original module weighting'),('drop_osmotic_desiccation_salt','Omit solute/ion module'),('equal_osmotic_submodules','Equal solute/ion subgroup weights')]):
 z=(eq if i==2 else d);z=z[(z.branch=='core65') & ~z.strict & (z.scheme=='lineage_BB') & (z.comparison==comp)].iloc[0];ax[2].plot([z.delta_q025*100,z.delta_q975*100],[i]*2,c='#1b6a82',lw=2);ax[2].plot(z.delta_median*100,i,'o',c='#1b6a82');ax[2].text(z.delta_q025*100,i-.14,label,fontsize=8)
ax[2].set(yticks=[],xlabel='Maintenance − energy partial R² (pp)',title='C  Panel sensitivity (core65)');ax[2].set_ylim(2.6,-.6);ax[2].axvline(0,c='.6',ls='--',lw=.8)
fig.text(.01,.025,'Intervals: 95% reweighting ranges; points: medians. Conditional on observed catalogues and coarse lineage blocks.\nEqual subgroup weighting is a post-result diagnostic, not a replacement primary endpoint.',fontsize=8)
fig.tight_layout(rect=[0,.105,1,1],w_pad=2)
for ext in ['png','pdf']:fig.savefig(H/f'portfolio_uncertainty_diagnostic.{ext}',dpi=200)
