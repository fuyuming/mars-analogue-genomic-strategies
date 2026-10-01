#!/usr/bin/env python3
"""Plot frozen joint-state diagnostics; expectations are descriptive, not fitted effects."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--input-dir', type=Path, required=True)
    p.add_argument('--output-prefix', type=Path, required=True)
    a = p.parse_args()
    s = pd.read_csv(a.input_dir/'joint_summary.tsv', sep='\t')
    g = pd.read_csv(a.input_dir/'joint_support.tsv', sep='\t')
    cohorts = ['original_primary','original_strict','new_desert_primary','new_desert_strict']
    labels = ['Original primary','Original strict','Dryland primary','Dryland strict']
    systems = ['KdpABC','ProVWX','EctABC','BetAB']
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5), gridspec_kw={'width_ratios':[1,1,1.1]})
    colors = ['#343B45','#DB9C44','#428E9D']
    for k,(cohort,label) in enumerate(zip(cohorts,labels)):
        ax=axes[k//2,k%2]
        d=s.set_index(['cohort','pair']).loc[cohort].reindex([v+'__CoxCut' for v in systems])
        x=np.arange(4)
        fields=['both_family_assigned','expected_both_pooled_family_assigned','expected_both_family_conditioned']
        for j,(field,legend) in enumerate(zip(fields,['Observed','Pooled expectation','Within-family expectation'])):
            ax.bar(x+(j-1)*.25,d[field],width=.23,color=colors[j],label=legend)
        ax.set_xticks(x,systems)
        ax.set_ylabel('MAGs with both marker sets')
        ax.set_title(f'{chr(65+k)}  {label} (family assigned: {int(d.family_assigned_n.iloc[0]):,})',loc='left',fontsize=10)
        ax.spines[['right','top']].set_visible(False)
        if k==0: ax.legend(frameon=False,fontsize=8)
    ax=axes[0,2]
    mat=np.array([[int(g[(g.cohort==c)&(g.pair==v+'__CoxCut')&(g.min_per_joint_state==3)].eligible_families.iloc[0]) for v in systems] for c in cohorts])
    ax.imshow(mat,cmap='Blues',vmin=0,vmax=6,aspect='auto')
    for (i,j),v in np.ndenumerate(mat):ax.text(j,i,str(v),ha='center',va='center',color='white' if v>=4 else '#222222',fontsize=14)
    ax.set_xticks(range(4),systems,rotation=25,ha='right');ax.set_yticks(range(4),labels,fontsize=9)
    ax.set_title('E  Families passing joint-state support',loc='left',fontsize=10)
    ax=axes[1,2];ax.axis('off')
    ax.text(0,1,'Support diagnostic',weight='bold',fontsize=12,va='top')
    ax.text(0,.87,'Each family must contain all four joint states\nwith ≥3 MAGs per state, plus ≥2 source/site\ngroups with ≥3 MAGs each.\n\nStrict dryland support remains zero at\n≥2, ≥3 and ≥5 MAGs per joint state.\n\nAll bars use the same family-assigned\nMAGs within each cohort. Expectations\nretain pooled or family-specific margins.\n\nCox/Cut: K03518–K03520 co-detection.\nNot a validated CO-oxidation endpoint.\nNot-all-detected is not confirmed absence.',fontsize=10,va='top',linespacing=1.5)
    fig.suptitle('Maintenance–Cox/Cut joint configurations: lineage composition and available support',fontsize=14,y=.99)
    fig.tight_layout(rect=[0,0,1,.96],w_pad=2.2,h_pad=3)
    a.output_prefix.parent.mkdir(parents=True,exist_ok=True)
    for ext in ['png','pdf','svg']:fig.savefig(str(a.output_prefix)+'.'+ext,dpi=220,bbox_inches='tight')
    plt.close(fig)

if __name__=='__main__':main()
