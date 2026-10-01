"""Produce complete posterior tables and descriptive diagnostic figures."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from run_bayesian import TRAITS

H=Path(__file__).resolve().parent;R=H/'results';F=H/'figures';F.mkdir(exist_ok=True)
GLOBAL=['ordinal_main','ordinal_shared','ordinal_prior_half','ordinal_prior_double','ordinal_baseline_diffuse_long']
ENV=['environment_main_long','environment_strict','environment_prior_half','environment_prior_double','environment_unpooled']
LABELS=['Surface / EPS','Osmotic / desiccation','Oxidative / redox','Trace-gas markers','Dormancy / resuscitation',
 'DNA repair','Cold / protein quality','Sulfur markers','Amino-acid routes','Cofactor routes']
SOURCES=['alaska_permafrost_reference','atacama_halite_rainfall_prjna484015','atacama_salt_crust_prjna351262',
 'australian_basalt_lava_tubes_bay_2025','mauna_loa_lava_tube_fishman_2023','qaidam_basin_mmag_1773','stordalen_mire_2019_hybrid_mags']
SL=['Alaska','Atacama halite','Atacama crust','Australia lava','Mauna Loa','Qaidam','Stordalen']
def read(tag,name):return pd.read_csv(R/tag/name,sep='\t')
def savefig(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(F/f'{name}.{ext}',dpi=200,bbox_inches='tight')
 plt.close(fig)
def md_table(df):
 cols=list(df.columns)
 return '\n'.join(['| '+' | '.join(cols)+' |','|'+'|'.join(['---']*len(cols))+'|']+[
  '| '+' | '.join(str(x).replace('|','/') for x in row)+' |' for row in df.itertuples(index=False,name=None)])

def main():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
 profiles=pd.concat([read(t,'source_profile_posteriors.tsv').assign(model=t) for t in GLOBAL],ignore_index=True)
 env=pd.concat([read(t,'environment_posteriors.tsv').assign(model=t) for t in ENV],ignore_index=True)
 profiles.to_csv(R/'all_source_posteriors.tsv',sep='\t',index=False);env.to_csv(R/'all_environment_posteriors.tsv',sep='\t',index=False)
 diags=[];pprows=[]
 for tag in GLOBAL+ENV+['simulation','ordinal_simulation','global_main','environment_main','ordinal_baseline_diffuse']:
  diag=json.loads((R/tag/'diagnostics.json').read_text());diags.append(dict(model=tag,**{k:v for k,v in diag.items() if k!='stat_variables'}))
  pp=read(tag,'posterior_predictive_checks.tsv');allpp=pp[pp.group.eq('ALL')]
  row=dict(model=tag,overall_outside95=int((~allpp.observed_in_95).sum()),overall_checks=len(allpp),
   all_groups_outside95=int((~pp.observed_in_95).sum()),all_group_checks=len(pp))
  cat=R/tag/'posterior_predictive_categories.tsv'
  if cat.exists():
   c=pd.read_csv(cat,sep='\t');ca=c[c.group.eq('ALL')];row.update(category_all_outside95=int((~ca.observed_in_95).sum()),category_all_checks=len(ca),category_group_outside95=int((~c.observed_in_95).sum()),category_group_checks=len(c))
  pprows.append(row)
 pd.DataFrame(diags).to_csv(R/'diagnostics_summary.tsv',sep='\t',index=False)
 pd.DataFrame(pprows).to_csv(R/'predictive_check_summary.tsv',sep='\t',index=False)
 # Complete main-cohort tables: all endpoints, no probability-based row filtering.
 blocks=['# 全终点后验表','', '差值单位：百分点；区间为边际95%后验分位区间。概率不作跨终点联合事件解释。参考均值以零协变量/零科或地点效应计算。','', '## 全域来源相对参考差值']
 a=profiles[profiles.model.eq('ordinal_main')]
 table=pd.DataFrame({'来源':a.source,'终点':a.trait,'中位差':(100*a['median']).round(2),
  '95%下限':(100*a.low95).round(2),'95%上限':(100*a.high95).round(2),'P正向':a.p_positive.round(3),
  'P差>5pp':a.p_gt_5pp.round(3),'P差<−5pp':a.p_lt_minus_5pp.round(3)})
 blocks.append(md_table(table))
 for tag in ['environment_main_long','environment_strict','environment_unpooled']:
  a=env[env.model.eq(tag)];blocks+=['',f'## Qaidam每倍增：{tag}','']
  table=pd.DataFrame({'压力/层级':a.exposure,'终点':a.trait,'中位差':(100*a['median']).round(2),
   '95%下限':(100*a.low95).round(2),'95%上限':(100*a.high95).round(2),'P正向':a.p_positive.round(3),'P绝对差<5pp':a.p_abs_lt_5pp.round(3)})
  blocks.append(md_table(table))
 (H/'complete_posterior_tables_zh.md').write_text('\n'.join(blocks)+'\n')
 # Heatmaps: all traits, sample sizes and uncertainty indicated without binary decisions.
 fig,axs=plt.subplots(1,2,figsize=(13,6),layout='constrained')
 limit=max(abs(profiles['median']).max()*100,1)
 for ax,tag,title in zip(axs,GLOBAL[:2],['All strict representatives','Families found in ≥2 sources']):
  a=profiles[profiles.model.eq(tag)];mat=a.pivot(index='trait',columns='source',values='median').loc[TRAITS,SOURCES]*100
  im=ax.imshow(mat,vmin=-limit,vmax=limit,cmap='RdBu_r',aspect='auto')
  for j,s in enumerate(SOURCES):
   for i,t in enumerate(TRAITS):
    r=a[a.source.eq(s)&a.trait.eq(t)].iloc[0];mark='•' if r.low95>0 or r.high95<0 else ''
    ax.text(j,i,f'{100*r["median"]:.1f}{mark}',ha='center',va='center',fontsize=8,color='white' if abs(100*r['median'])>.65*limit else 'black')
  counts=a.groupby('source').n_source.first();ax.set_xticks(range(7),[f'{l}\nn={counts[s]}' for l,s in zip(SL,SOURCES)],rotation=45,ha='right')
  ax.set_yticks(range(10),LABELS);ax.set_title(title,loc='left',fontweight='bold')
 fig.colorbar(im,ax=axs,shrink=.6,label='Reference mean support shift (percentage points)')
 fig.suptitle('Source-associated genomic profiles | ordinal model\n• Marginal 95% interval excludes zero; source includes study background',fontsize=11)
 savefig(fig,'source_profiles')
 # Within-site slopes: all three exposures and all ten traits, three model assumptions.
 fig,axs=plt.subplots(1,3,figsize=(13,6),sharey=True,layout='constrained')
 styles=[('environment_main_long','Hierarchical main','#176d9c',-.18),('environment_strict','Strict-quality subset','#bd5c27',0),('environment_unpooled','No cross-trait shrinkage','#666666',.18)]
 for ax,e in zip(axs,['water','EC','TOC']):
  for tag,label,color,off in styles:
   a=env[env.model.eq(tag)&env.exposure.eq('log2_'+e+'_within')].set_index('trait').loc[TRAITS]
   x=a['median'].to_numpy()*100;lo=a.low95.to_numpy()*100;hi=a.high95.to_numpy()*100
   ax.errorbar(x,np.arange(10)+off,xerr=[x-lo,hi-x],fmt='o',ms=3,lw=1,color=color,label=label)
  ax.axvline(0,c='.7',lw=.8);ax.set_title(e+' within site',loc='left',fontweight='bold');ax.set_xlabel('Change per doubling (pp)');ax.grid(axis='x',alpha=.15)
 selected=env[env.model.isin([s[0] for s in styles])&env.exposure.str.endswith('_within')]
 limit=100*max(selected.low95.abs().max(),selected.high95.abs().max())*1.08
 for ax in axs:ax.set_xlim(-limit,limit)
 axs[0].set_yticks(range(10),LABELS);axs[0].invert_yaxis()
 handles,labels=axs[-1].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.55,-.045),ncol=3,fontsize=8)
 fig.suptitle('Qaidam pressure associations | median and marginal 95% interval\nRecovered-MAG support, not measured activity or a causal effect',fontsize=11)
 savefig(fig,'qaidam_pressure_associations')
 # Full count-category posterior checks at global level.
 c=read('ordinal_main','posterior_predictive_categories.tsv');fig,axs=plt.subplots(2,5,figsize=(14,5.5),layout='constrained')
 for ax,t,label in zip(axs.flat,TRAITS,LABELS):
  a=c[c.group.eq('ALL')&c.trait.eq(t)].sort_values('category');x=a.category.to_numpy()
  ax.bar(x,a.observed,width=.7,color='#284759',alpha=.8,label='Observed')
  ax.errorbar(x,a.predictive_median,yerr=[a.predictive_median-a.predictive_low95,a.predictive_high95-a.predictive_median],fmt='o',ms=2,color='#dd8e26',lw=.9,label='Posterior predictive 95%')
  ax.set_title(label,fontsize=9,loc='left');ax.set_xlabel('Support count');ax.set_ylabel('Fraction');ax.set_ylim(bottom=0);ax.xaxis.set_major_locator(MaxNLocator(integer=True,nbins=5))
 axs[0,0].legend(fontsize=6);fig.suptitle('Ordinal count-distribution checks | all 879 genomes',fontsize=12)
 savefig(fig,'ordinal_count_checks')
 # Predetermined trace coordinates, not only best-looking chains.
 for tag in ['ordinal_main','environment_main_long']:
  tr=np.load(R/tag/'trace_core.npz');fig,axs=plt.subplots(2,3,figsize=(12,4),layout='constrained')
  for col,j in enumerate([0,4,9]):
   for chain in range(4):
    axs[0,col].plot(tr['beta'][chain,:,0,j],lw=.4,alpha=.65)
   axs[0,col].set_title('beta[0]: '+LABELS[j],fontsize=8)
   v=tr['scale'];k=j if v.shape[-1]==10 else col
   for chain in range(4):axs[1,col].plot(v[chain,:,k],lw=.4,alpha=.65)
   axs[1,col].set_title(f'scale[{k}]',fontsize=8);axs[1,col].set_xlabel('Posterior draw')
  fig.suptitle(tag+' | four-chain traces');savefig(fig,'trace_'+tag)
 print(pd.DataFrame(diags)[['model','rhat_max_core','ess_bulk_min_core','ess_tail_min_core','divergences','core_gate_pass']].to_string(index=False))
 print(pd.DataFrame(pprows).to_string(index=False))

if __name__=='__main__':main()
