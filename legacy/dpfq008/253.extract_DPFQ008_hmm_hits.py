#!/usr/bin/env python3
"""Extract catalogue-wide DPFQ008 HMM hits and audit CXXCH motifs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


MOTIF_RE = re.compile(r"C..CH")
THRESHOLDS = {
    "wide": (1e-3, 0.50),
    "moderate": (1e-5, 0.60),
    "strict": (1e-10, 0.70),
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hits", required=True, type=Path)
    p.add_argument("--catalogue-fasta", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args()


def best_hits(path: Path) -> dict[str, dict[str, str]]:
    best: dict[str, dict[str, str]] = {}
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            old = best.get(row["target"])
            if old is None or float(row["domain_ievalue"]) < float(old["domain_ievalue"]):
                best[row["target"]] = row
    return best


def extract(path: Path, wanted: set[str]) -> dict[str, str]:
    found: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []

    def commit() -> None:
        if current is None:
            return
        seq = "".join(chunks).rstrip("*")
        if not seq:
            raise RuntimeError(f"Empty sequence for selected hit: {current}")
        if current in found and found[current] != seq:
            raise RuntimeError(f"Conflicting duplicate: {current}")
        found[current] = seq

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                commit()
                identifier = line[1:].split()[0]
                current = identifier if identifier in wanted else None
                chunks = []
            elif current is not None:
                chunks.append(line.strip())
        commit()
    missing = wanted.difference(found)
    if missing:
        raise RuntimeError(f"Missing selected sequences: {len(missing)}; first={sorted(missing)[:3]}")
    return found


def main() -> int:
    a = parse_args()
    a.hits.resolve(strict=True)
    a.catalogue_fasta.resolve(strict=True)
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    best = best_hits(a.hits)
    wide = {
        protein for protein, row in best.items()
        if float(row["domain_ievalue"]) <= THRESHOLDS["wide"][0]
        and float(row["target_coverage"]) >= THRESHOLDS["wide"][1]
    }
    if len(wide) != 126:
        raise RuntimeError(f"Expected frozen wide set of 126, saw {len(wide)}")
    sequences = extract(a.catalogue_fasta, wide)

    a.outdir.mkdir(parents=True)
    fasta = a.outdir / "DPFQ008_HMM_wide_126.faa"
    with fasta.open("x", encoding="utf-8") as handle:
        for protein in sorted(sequences):
            handle.write(f">{protein}\n")
            seq = sequences[protein]
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i+80] + "\n")

    fields = [
        "protein_id", "dataset_id", "genome_id", "original_protein_id",
        "length_aa", "domain_ievalue", "full_evalue", "target_coverage",
        "hmm_coverage", "wide", "moderate", "strict", "cxxch_count",
        "cxxch_positions", "cxxch_motifs",
    ]
    rows = []
    for protein in sorted(wide):
        row = best[protein]
        seq = sequences[protein]
        motifs = list(MOTIF_RE.finditer(seq))
        dataset, genome, original = protein.split("|", 2)
        membership = {
            name: float(row["domain_ievalue"]) <= evalue and float(row["target_coverage"]) >= coverage
            for name, (evalue, coverage) in THRESHOLDS.items()
        }
        rows.append({
            "protein_id": protein, "dataset_id": dataset, "genome_id": genome,
            "original_protein_id": original, "length_aa": len(seq),
            "domain_ievalue": row["domain_ievalue"], "full_evalue": row["full_evalue"],
            "target_coverage": row["target_coverage"], "hmm_coverage": row["hmm_coverage"],
            **{name: int(value) for name, value in membership.items()},
            "cxxch_count": len(motifs),
            "cxxch_positions": ";".join(str(m.start()+1) for m in motifs),
            "cxxch_motifs": ";".join(m.group(0) for m in motifs),
        })
    with (a.outdir / "DPFQ008_HMM_hit_architecture.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    lines = ["DPFQ008 HMM hit extraction and motif audit"]
    for name in THRESHOLDS:
        selected = [row for row in rows if row[name]]
        lines.extend([
            f"{name}_proteins={len(selected)}",
            f"{name}_datasets={len({row['dataset_id'] for row in selected})}",
            f"{name}_with_2plus_CXXCH={sum(row['cxxch_count'] >= 2 for row in selected)}",
            f"{name}_with_any_CXXCH={sum(row['cxxch_count'] >= 1 for row in selected)}",
        ])
    lines.extend([
        f"fasta_sha256={hashlib.sha256(fasta.read_bytes()).hexdigest()}",
        "claim_boundary=custom HMM plus motif conservation supports family distribution, not biochemical activity",
        "status=PASS",
    ])
    (a.outdir / "run_summary.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
