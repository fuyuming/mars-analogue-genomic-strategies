from pathlib import Path
import copy,json,hashlib
import pandas as pd,numpy as np
from Bio import Phylo
H=Path(__file__).resolve().parent;O=H/'host';meta=pd.read_csv(O/'strict1243_integrated_metadata.tsv',sep='\t').set_index('tip_id');rows=[];summaries=[];sens=[]
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
for marker,domain,pref in [('bac120','Bacteria','B'),('ar53','Archaea','A')]:
 tree=Phylo.read(O/f'tree_completed/inference/strict/{marker}_LG_F_R10.treefile','newick');df=meta.loc[[t.name for t in tree.get_terminals()]]
 # Both roots are visual/model conventions, not inferred ancestors.
 alternate=copy.deepcopy(tree);alternate.root_with_outgroup(min(df.index));ab=partition(alternate,df)
 tree.root_at_midpoint();tree.ladderize();blocks=partition(tree,df,unrooted=True);assert {frozenset(x[3]) for x in blocks}=={frozenset(x[3]) for x in partition(alternate,df,unrooted=True)};bs={frozenset(x[3]) for x in blocks};als={frozenset(x[3]) for x in ab};sens.append(dict(domain=domain,unrooted_edge_blocks=len(bs),lex_tip_root_blocks=len(als),identical_member_blocks=len(bs&als),genomes_in_identical_blocks=sum(len(x) for x in bs&als),note='Unrooted edge-defined partition verified identical under both roots; rooted-only lex-tip count is diagnostic'))
 keep=[];rename={}
 for i,(path,rank,name,ids) in enumerate(blocks,1):
  lid=f'{pref}63_L{i:03d}';rep=min(ids);keep.append(rep);rename[rep]=lid
  for tip in ids:rows.append(dict(tip_id=tip,genome_id=df.loc[tip,'genome_id'],domain=domain,lineage_id=lid,collapsed_rank=rank,collapsed_name=name,lineage_path=path))
  summaries.append(dict(domain=domain,lineage_id=lid,collapsed_rank=rank,collapsed_name=name,genomes=len(ids),sources=df.loc[ids,'dataset_id'].nunique(),display_representative=rep))
 collapsed=copy.deepcopy(tree)
 for tip in list(collapsed.get_terminals()):
  if tip.name not in keep:collapsed.prune(tip)
 for tip in collapsed.get_terminals():tip.name=rename[tip.name]
 for cl in collapsed.get_nonterminals():cl.name=None;cl.confidence=None
 Phylo.write(collapsed,O/f'{marker}_collapsed_display.newick','newick');Phylo.write(tree,O/f'{marker}_midpoint_display.newick','newick')
m=pd.DataFrame(rows);assert len(m)==1243 and m.tip_id.is_unique;m.to_csv(O/'strict_tree_lineage_mapping.tsv',sep='\t',index=False);pd.DataFrame(summaries).to_csv(O/'strict_lineage_summary.tsv',sep='\t',index=False);pd.DataFrame(sens).to_csv(O/'root_partition_sensitivity.tsv',sep='\t',index=False);print(pd.DataFrame(sens).to_string(index=False));print(pd.DataFrame(summaries).groupby(['domain','collapsed_rank']).size())
