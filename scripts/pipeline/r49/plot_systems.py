from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np,pandas as pd
H=Path(__file__).resolve().parent/'host';s=pd.read_csv(H/'submodule_decomposition.tsv',sep='\t');z=pd.read_csv(H/'system_state_summary.tsv',sep='\t');g=pd.read_csv(H/'lineage_group_support.tsv',sep='\t');plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,3,figsize=(13.6,4.9),gridspec_kw={'width_ratios':[1.08,1.25,1]})
a=s[(s.branch=='core65') & ~s.strict];labs={'betaine':'Betaine marker','compatible_solute_transport':'Compatible-solute transport','ectoine':'Ectoine markers','potassium_homeostasis':'Kdp pump + regulator','trehalose':'Trehalose-related markers'}
ax[0].barh(range(len(a)),100*a.source_partial_r2,color=['#8ca9b5' if k!='potassium_homeostasis' else '#1b6a82' for k in a.submodule]);ax[0].set(yticks=range(len(a)),yticklabels=[labs[k] for k in a.submodule],xlabel='Source partial R² (%)',title='A  Original subgroup structure');ax[0].invert_yaxis()
for i,sys in enumerate(['KdpABC','ProVWX','EctABC','BetAB']):
 for cohort,off,col in [('original_strict',-.14,'#1b6a82'),('new_desert_strict',.14,'#b26b36')]:
  r=z[(z.cohort==cohort)&(z.system==sys)].iloc[0];ax[1].plot([r.conditional_reweight_q025,r.conditional_reweight_q975],[i+off]*2,c=col,lw=1.8);ax[1].plot(r.complete_given_any,i+off,'o',c=col,ms=4,label=('Original: lineage weights' if off<0 else 'Desert: site weights') if i==0 else None);ax[1].text(1.04,i+off,f'{r.panel_complete}/{r.panel_complete+r.partial}',va='center',color=col,fontsize=8)
ax[1].set(yticks=range(4),yticklabels=['KdpABC core','ProVWX components','EctABC core','BetAB pair'],xlabel='All panel markers / any panel marker',title='B  Co-detection in strict-quality MAGs',xlim=(-.03,1.26),xticks=[0,.25,.5,.75,1],ylim=(4.5,-.6));ax[1].legend(frameon=False,fontsize=8,loc='lower left')
for i,rank in enumerate(['class','order','family']):
 for cohort,off,col in [('new_desert_primary',-.18,'#1b6a82'),('new_desert_strict',.18,'#b26b36')]:
  r=g[(g.cohort==cohort)&(g['rank']==rank)].iloc[0];ax[2].barh(i+off,r.n,height=.3,color=col,label=('Primary' if off<0 else 'Strict') if i==0 else None);ax[2].text(r.n+6,i+off,str(r.n),va='center',fontsize=8)
ax[2].set(yticks=range(3),yticklabels=['Class','Order','Family'],xlabel='MAGs in supported lineages',title='C  Desert within-lineage support',xlim=(0,650));ax[2].invert_yaxis();ax[2].legend(frameon=False,fontsize=8,loc='lower right')
fig.text(.012,.032,'B: 95% conditional reweighting ranges; denominators are MAGs with at least one marker. Co-detection is not pathway activity.\nC: At least two sites with ≥5 MAG representatives each in a lineage; genomic units are not independent environmental replicates.',fontsize=8)
fig.tight_layout(rect=[0,.11,1,1],w_pad=2.4)
for ext in ['png','pdf']:fig.savefig(H/f'system_configuration_diagnostic.{ext}',dpi=200)
