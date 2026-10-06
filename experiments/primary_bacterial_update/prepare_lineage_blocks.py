from pathlib import Path
import copy
import pandas as pd
from Bio import Phylo
H=Path(__file__).resolve().parent;O=H/'host'
meta=pd.read_csv(O/'primary4465_integrated_metadata.tsv',sep='\t').set_index('tip_id')
df=meta[meta.expected_domain.eq('Bacteria')];tree=Phylo.read(O/'primary_bac120.treefile','newick')
assert len(tree.get_terminals())==4088 and set(t.name for t in tree.get_terminals())==set(df.index)
ranks=['phylum','class','order','family','genus']
def partition(tree,df,unrooted=False):
 desc={frozenset(t.name for t in c.get_terminals()):c for c in tree.find_clades()};out=[]
 if unrooted:
  allids=frozenset(df.index);desc.update({allids-key:cl for key,cl in list(desc.items()) if allids-key})
 def split(ids,rankidx,path):
  rank=ranks[rankidx]
  for value,g in df.loc[ids].groupby(rank,sort=True):
   names=list(g.index);key=frozenset(names);label=path+f'/{rank}:{value}'
   if key in desc:out.append((label,rank,value,names))
   elif rankidx<4:split(names,rankidx+1,label)
   else:
    for name in names:out.append((label+'/genome:'+name,'genome',name,[name]))
 split(list(df.index),0,'');return out
blocks=partition(tree,df,unrooted=True);alternate=copy.deepcopy(tree);alternate.root_with_outgroup(min(df.index))
assert {frozenset(v[3]) for v in blocks}=={frozenset(v[3]) for v in partition(alternate,df,unrooted=True)}
rows=[]
for k,(path,rank,name,ids) in enumerate(blocks,1):
 for tip in ids:rows.append(dict(tip_id=tip,genome_id=df.loc[tip,'genome_id'],domain='Bacteria',lineage_id=f'B69_L{k:03d}',collapsed_rank=rank,collapsed_name=name,lineage_path=path))
pd.DataFrame(rows).to_csv(O/'primary_tree_lineage_mapping.tsv',sep='\t',index=False)
print('Primary bacterial unrooted groups:',len(blocks))
