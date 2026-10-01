"""Read existing remote mapper and GFF files; stream only predeclared targets to stdout.

The annotation/prodigal root is supplied via --dbroot (no private absolute path
is embedded). The extraction algorithm is unchanged from the authoring version.
"""
from pathlib import Path
from collections import defaultdict
import argparse,csv,sys,re
_ap=argparse.ArgumentParser(description="Stream predeclared dryland target proteins (Kofam mapper + Prodigal GFF)")
_ap.add_argument('--dbroot',type=Path,
                 default=Path(__import__('os').environ.get('NEE_DESERT_DBROOT','desert_annotation_root')),
                 help="directory holding full_kofam_1.3.0/ and prodigal/ (or set NEE_DESERT_DBROOT)")
R=_ap.parse_args().dbroot
kos=set('K01546 K01547 K01548 K07646 K07667 K02000 K02001 K02002 K06718 K00836 K06720 K10674 K00108 K00130'.split())
hits=defaultdict(dict)
with (R/'full_kofam_1.3.0/mapper_one_line.tsv').open() as f:
 for line in f:
  a=line.rstrip().split('\t');selected=kos.intersection(a[1:])
  if selected:
   genome,protein=a[0].split('|',1);hits[genome][protein]=selected
fields=['genome_id','accession','protein_id','KO','contig','cds_rank','start','end','strand','partial','contig_length']
w=csv.DictWriter(sys.stdout,fieldnames=fields,delimiter='\t');w.writeheader()
for genome in sorted(hits):
 found=set();length=0
 with (R/'prodigal'/genome/'genes.gff').open() as f:
  for line in f:
   if line.startswith('# Sequence Data:'):length=int(re.search(r'seqlen=(\d+)',line).group(1))
   if line.startswith('#'):continue
   a=line.rstrip().split('\t')
   if len(a)!=9 or a[2]!='CDS':continue
   attr=dict(x.split('=',1) for x in a[8].split(';') if '=' in x)
   rank=int(attr['ID'].rsplit('_',1)[1]);pid=a[0]+'_'+str(rank)
   if pid not in hits[genome]:continue
   found.add(pid)
   for ko in sorted(hits[genome][pid]):
    w.writerow(dict(genome_id=genome,accession=genome,protein_id=pid,KO=ko,contig=a[0],cds_rank=rank,start=a[3],end=a[4],strand=a[6],partial=attr['partial'],contig_length=length))
 assert found==set(hits[genome]),(genome,'unmapped target protein')
