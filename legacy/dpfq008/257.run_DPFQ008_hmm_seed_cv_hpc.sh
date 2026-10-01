#!/usr/bin/env bash
set -euo pipefail

family_faa=${1:?20-member family FASTA required}
catalogue_faa=${2:?catalogue FASTA required}
full_hit_table=${3:?full-profile hit architecture table required}
outdir=${4:?output directory required}
threads=${5:-32}

[[ -s "$family_faa" && -s "$catalogue_faa" && -s "$full_hit_table" ]] || { echo "Input absent" >&2; exit 2; }
[[ ! -e "$outdir" ]] || { echo "Refusing to overwrite: $outdir" >&2; exit 3; }
tmp="${outdir}.tmp.$$"
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT

python3 - "$family_faa" "$tmp" <<'PY'
import sys
from pathlib import Path
source, out = Path(sys.argv[1]), Path(sys.argv[2])
records=[]
head=None; seq=[]
def commit():
    if head is not None: records.append((head, ''.join(seq)))
with source.open() as h:
    for line in h:
        if line.startswith('>'):
            commit(); head=line[1:].rstrip('\n'); seq=[]
        else: seq.append(line.strip())
    commit()
if len(records)!=20 or any(not s for _,s in records): raise SystemExit(f'family gate failed: {len(records)}')
for fold in range(1,6):
    held=[r for i,r in enumerate(records,1) if (i-1)%5+1==fold]
    train=[r for i,r in enumerate(records,1) if (i-1)%5+1!=fold]
    if len(train)!=16 or len(held)!=4: raise SystemExit('fold-size gate failed')
    with (out/f'fold{fold}.train.faa').open('x') as h:
        for head,seq in train: h.write(f'>{head}\n{seq}\n')
    with (out/f'fold{fold}.heldout_ids.txt').open('x') as h:
        for head,_ in held:
            fields=head.split(maxsplit=1)
            if len(fields)!=2 or '|' not in fields[1]: raise SystemExit(f'header mapping gate failed: {head}')
            h.write(fields[1]+'\n')
PY

for fold in 1 2 3 4 5; do
  mafft --auto --thread "$threads" "$tmp/fold${fold}.train.faa" \
    > "$tmp/fold${fold}.mafft.faa" 2> "$tmp/fold${fold}.mafft.log"
  hmmbuild "$tmp/fold${fold}.hmm" "$tmp/fold${fold}.mafft.faa" \
    > "$tmp/fold${fold}.hmmbuild.log"
  hmmsearch --cpu "$threads" --noali --domtblout "$tmp/fold${fold}.domtblout" \
    -E 1e-3 --domE 1e-3 "$tmp/fold${fold}.hmm" "$catalogue_faa" \
    > "$tmp/fold${fold}.hmmsearch.txt"
done

python3 - "$full_hit_table" "$tmp" <<'PY'
import csv, sys
from pathlib import Path
full_table, root = Path(sys.argv[1]), Path(sys.argv[2])
with full_table.open() as h:
    full={r['protein_id'] for r in csv.DictReader(h,delimiter='\t') if r['strict']=='1'}
if len(full)!=63: raise SystemExit(f'full strict gate failed: {len(full)}')
fold_sets=[]; summaries=[]
for fold in range(1,6):
    best={}
    with (root/f'fold{fold}.domtblout').open() as h:
        for line in h:
            if line.startswith('#'): continue
            f=line.split()
            target=f[0]; tlen=int(f[2]); ie=float(f[12]); cov=(int(f[18])-int(f[17])+1)/tlen
            old=best.get(target)
            if old is None or ie<old[0]: best[target]=(ie,cov)
    strict={p for p,(ie,cov) in best.items() if ie<=1e-10 and cov>=.70}
    held={x.strip() for x in (root/f'fold{fold}.heldout_ids.txt').read_text().splitlines() if x.strip()}
    fold_sets.append(strict)
    summaries.append({
        'fold':fold,'training_seeds':16,'heldout_seeds':4,
        'heldout_recovered':len(held & strict),'strict_hits':len(strict),
        'overlap_full63':len(strict & full),'union_full63':len(strict | full),
        'jaccard_vs_full63':len(strict & full)/len(strict | full),
        'full63_recovered':len(strict & full)/63,
    })
with (root/'fold_summary.tsv').open('x',newline='') as h:
    w=csv.DictWriter(h,fieldnames=list(summaries[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(summaries)
all_hits=sorted(set().union(*fold_sets))
with (root/'hit_fold_frequency.tsv').open('x',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['protein_id','folds_detected','in_full_strict63'])
    for p in all_hits: w.writerow([p,sum(p in s for s in fold_sets),int(p in full)])
consensus={p for p in all_hits if sum(p in s for s in fold_sets)>=4}
lines=[
    'DPFQ008 five-fold seed-profile stability',
    'folds=5','training_seeds_per_fold=16','heldout_seeds_per_fold=4',
    f'heldout_recovered_total={sum(x["heldout_recovered"] for x in summaries)}/20',
    f'strict_hits_min={min(x["strict_hits"] for x in summaries)}',
    f'strict_hits_max={max(x["strict_hits"] for x in summaries)}',
    f'jaccard_vs_full63_min={min(x["jaccard_vs_full63"] for x in summaries):.6f}',
    f'jaccard_vs_full63_max={max(x["jaccard_vs_full63"] for x in summaries):.6f}',
    f'hits_in_4plus_folds={len(consensus)}',
    f'full63_in_4plus_folds={len(consensus & full)}/63',
    'claim_boundary=cross-seed stability bounds operational-profile dependence but is not an external specificity benchmark',
    'status=PASS'
]
(root/'run_summary.txt').write_text('\n'.join(lines)+'\n')
print('\n'.join(lines))
PY

mv "$tmp" "$outdir"
trap - EXIT
cat "$outdir/fold_summary.tsv"
