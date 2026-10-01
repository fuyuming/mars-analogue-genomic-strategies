from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np,csv
H=Path(__file__).resolve().parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
def save(fig,name):
 fig.savefig(H/'assets'/f'{name}.png',dpi=220,bbox_inches='tight',facecolor='white')
 fig.savefig(H/'assets'/f'{name}.pdf',bbox_inches='tight',facecolor='white');plt.close(fig)
fig,ax=plt.subplots(figsize=(10,5.7));ax.set(xlim=(0,10),ylim=(0,5.7));ax.axis('off')
def box(x,y,w,h,title,body,color):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.10',fc=color,ec='#9AABB5',lw=.8))
 ax.text(x+.13,y+h-.22,title,fontsize=11,weight='bold',va='top',color='#173647')
 ax.text(x+.13,y+h-.62,body,fontsize=9.5,va='top',linespacing=1.45,color='#263942')
box(.15,3.5,4.4,1.6,'CONSTRAINTS','Hydration • salinity / ion chemistry • pH\nTemperature • exposure / shielding\nPreserve measured units and sampling scale','#EAF1F7')
box(5.35,3.5,4.4,1.6,'RESOURCE CONTEXT','Organic-carbon stocks and composition\nEnergy donors / acceptors and mineral context\nStocks and elements do not establish usable flux','#F6F0DE')
box(.15,.98,9.6,1.4,'GENOMIC STRATEGY COMBINATIONS','Cellular maintenance  |  Energy acquisition  |  Carbon assimilation\nLineage and detection as competing explanations; strategies need not be mutually exclusive','#EAF3EE')
for x in [2.4,7.6]:ax.annotate('',xy=(x,2.48),xytext=(x,3.36),arrowprops={'arrowstyle':'->','color':'#46616D','lw':1.5,'linestyle':'--'})
ax.text(5,2.92,'Conditional relationships to test',ha='center',fontsize=10,color='#5B6266',bbox={'facecolor':'white','edgecolor':'none','pad':2})
ax.text(5,.42,'Bounded implications: geological targeting for life detection | candidate functions for engineered resource use',ha='center',fontsize=9)
ax.text(.1,5.5,'Conceptual framework — no causal effect or Mars performance estimated',fontsize=11,weight='bold')
save(fig,'figure1_framework')
rows=[['bacteria_core65',5.0064,2.1932],['bacteria_restored74',4.419,2.1932],['archaea_core65',8.2073,17.0083]]
# Restored maintenance 4.42 is reported rounded; do not imply unseen precision.
rows[1][1]=4.42
with (H/'figure2_source_partition.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['cohort','maintenance_partial_R2_percent','energy_partial_R2_percent']);w.writerows(rows)
fig,axes=plt.subplots(3,1,figsize=(8.5,9),gridspec_kw={'height_ratios':[1.3,1.2,1]});fig.subplots_adjust(hspace=.65,left=.29,right=.88,top=.96,bottom=.08)
a=axes[0];y=np.arange(3)
a.barh(y-.16,[r[1] for r in rows],height=.30,color='#286B84',label='Maintenance')
a.barh(y+.16,[r[2] for r in rows],height=.30,color='#CB9344',label='Energy')
a.set(yticks=y,yticklabels=['Bacteria · core65','Bacteria · restored74','Archaea · core65'],xlabel='Source partial R² (%)',xlim=(0,20));a.invert_yaxis();a.legend(frameon=False,loc='lower right');a.set_title('A   Source partitioning',loc='left',weight='bold')
for i,r in enumerate(rows):
 for j,v in enumerate(r[1:]):a.text(v+.15,i+(-.16 if j==0 else .16),f'{v:.2f}',va='center',fontsize=9)
b=axes[1];vals=[-.02813,-.02226,-.04130,-.03351];labels=['Core65 · primary','Restored74 · primary','Core65 · strict','Restored74 · strict']
b.scatter(vals,range(4),color='#286B84',s=42);b.axvline(0,color='0.65',lw=1);b.set(yticks=range(4),yticklabels=labels,xlim=(-.055,.008),xlabel='Energy − maintenance partial R²');b.invert_yaxis();b.set_title('B   Bacterial contrasts (point estimates)',loc='left',weight='bold')
for i,v in enumerate(vals):b.text(v-.0018,i,f'{v:.5f}',ha='right',va='center',fontsize=9)
c=axes[2];c.barh([0,1],[58/65*100,63/74*100],height=.48,color='#599B81');c.set(yticks=[0,1],yticklabels=['Core65','Restored74'],xlim=(0,100),xlabel='Targets detected in independent reads (%)');c.invert_yaxis();c.set_title('C   Occurrence evidence from one independent study',loc='left',weight='bold')
for i,(n,d) in enumerate([(58,65),(63,74)]):c.text(2,i,f'{n}/{d}',va='center',color='white',weight='bold')
save(fig,'figure2_existing_results')
vals=[('Main',17.70,14.12,21.16),('Shared families',17.04,12.66,21.22),('Prior half',17.14,13.61,20.65),('Prior double',17.93,14.33,21.42),('Diffuse baseline',17.43,14,20.77)]
with (H/'figure3_bayesian_source.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['model','median_pp','lower95_CrI_pp','upper95_CrI_pp']);w.writerows(vals)
fig,ax=plt.subplots(figsize=(8.5,4.2));fig.subplots_adjust(left=.25,right=.98,bottom=.2,top=.86)
for i,(lab,m,l,u) in enumerate(vals):
 ax.errorbar(m,i,xerr=[[m-l],[u-m]],fmt='o',color='#286B84',capsize=3);ax.text(22.1,i,f'{m:.2f} [{l:.2f}, {u:.2f}]',va='center',fontsize=9)
ax.axvline(0,c='0.65',lw=1);ax.set(yticks=range(5),yticklabels=[v[0] for v in vals],xlim=(-1,33),xlabel='Stordalen − Qaidam osmotic-marker support (percentage points)');ax.invert_yaxis();ax.set_title('Existing Bayesian source contrast · median and 95% credible interval',loc='left',weight='bold',fontsize=11)
save(fig,'figure3_bayesian_profiles')
