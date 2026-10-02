"""Bind frozen genome IDs to local plain FASTA paths without running analyses."""
from pathlib import Path
import argparse,csv,json,hashlib

def main():
 p=argparse.ArgumentParser();p.add_argument('--paths',type=Path,required=True,help='TSV: genome_id, fasta_path');p.add_argument('--config',type=Path,required=True);a=p.parse_args()
 cfg=json.loads(a.config.read_text());root=Path(cfg['root']);assert root.is_absolute() and not root.exists(),'Choose a fresh absolute output directory'
 src=Path(__file__).resolve().parent/'frozen_genome_manifest.tsv'
 rows=list(csv.DictReader(src.open(),delimiter='\t'));paths=list(csv.DictReader(a.paths.open(),delimiter='\t'))
 lookup={r['genome_id']:r['fasta_path'] for r in paths};assert len(lookup)==len(paths),'Duplicate genome mapping'
 assert set(lookup)=={r['genome_id'] for r in rows},'Missing or unexpected genomes'
 for row in rows:
  f=Path(lookup[row['genome_id']]).expanduser().resolve();assert f.is_file(),f
  with f.open('rb') as handle:assert handle.read(1)==b'>',f'Plain FASTA required: {f}'
  row['fasta_path']=str(f)
 assert len(rows)==cfg['expected_union']==4477
 root.mkdir(parents=True);dest=root/'union_genome_manifest.tsv'
 with dest.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 cfg['manifest_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();(root/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
 print(root/'config.json')

if __name__=='__main__':main()
