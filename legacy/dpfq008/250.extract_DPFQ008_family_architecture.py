#!/usr/bin/env python3
"""Extract the frozen DPFQ008 MMseqs50 family and audit haem-binding motifs.

This script is read-only with respect to the cluster table and source FASTA. It
refuses to overwrite an output directory and requires the frozen 20-member
family gate before writing results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from pathlib import Path


REPRESENTATIVE = (
    "qaidam_basin_mmag_1773|qaidam_basin_mmag_1773_CHGRbin21|"
    "k141_121601_length_11112_cov_23.0792_5"
)
EXPECTED_MEMBERS = 20


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-tsv", required=True, type=Path)
    parser.add_argument("--all-sequences", required=True, type=Path)
    parser.add_argument("--mag-manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    return parser.parse_args()


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def split_id(value: str) -> tuple[str, str, str]:
    fields = value.split("|", 2)
    if len(fields) != 3:
        raise ValueError(f"Malformed protein identifier: {value}")
    return fields[0], fields[1], fields[2]


def family_members(path: Path) -> list[str]:
    members: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed cluster row {line_number}")
            if fields[0] == REPRESENTATIVE:
                members.append(fields[1])
    if len(members) != EXPECTED_MEMBERS:
        raise RuntimeError(f"Expected {EXPECTED_MEMBERS} family members, observed {len(members)}")
    if len(set(members)) != EXPECTED_MEMBERS:
        raise RuntimeError("Duplicate DPFQ008 family members")
    return members


def read_fasta_subset(path: Path, wanted: set[str]) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    parts: list[str] = []

    def commit() -> None:
        if current in wanted:
            sequence = "".join(parts).upper()
            # The frozen concatenated FASTA contains an empty alias header
            # immediately before some metadata-bearing records. Ignore only
            # those zero-length aliases; a member must still resolve to one
            # non-empty sequence by the final gate below.
            if not sequence:
                return
            if re.search(r"[^ACDEFGHIKLMNPQRSTVWY]", sequence):
                raise ValueError(f"Non-standard sequence: {current}")
            if current in sequences and sequences[current] != sequence:
                raise ValueError(f"Conflicting duplicate FASTA records: {current}")
            sequences[current] = sequence

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                commit()
                current = line[1:].split()[0]
                parts = []
            elif current is not None:
                parts.append(line.strip())
        commit()
    missing = sorted(wanted - sequences.keys())
    if missing:
        raise RuntimeError("Family sequences absent from FASTA:\n" + "\n".join(missing))
    return sequences


def manifest_rows(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    result = {(row["dataset_id"], row["genome_id"]): row for row in rows}
    if len(result) != len(rows):
        raise RuntimeError("Duplicate dataset/genome rows in MAG manifest")
    return result


def main() -> int:
    args = arguments()
    for path in (args.cluster_tsv, args.all_sequences, args.mag_manifest):
        path.resolve(strict=True)
    if args.outdir.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.outdir}")

    members = family_members(args.cluster_tsv)
    sequences = read_fasta_subset(args.all_sequences, set(members))
    manifest = manifest_rows(args.mag_manifest)

    rows = []
    for index, protein_id in enumerate(sorted(members), 1):
        dataset, genome, original = split_id(protein_id)
        meta = manifest.get((dataset, genome))
        if meta is None:
            raise RuntimeError(f"Carrier absent from primary MAG manifest: {dataset}|{genome}")
        sequence = sequences[protein_id]
        motifs = [(match.start() + 1, match.group()) for match in re.finditer(r"C..CH", sequence)]
        rows.append({
            "stable_id": f"DPFQ008_{index:02d}",
            "protein_id": protein_id,
            "dataset_id": dataset,
            "genome_id": genome,
            "original_protein_id": original,
            "ani_cluster_95": meta["ani_cluster_95"],
            "gtdb_taxonomy_r226": meta["gtdb_taxonomy_r226"],
            "length_aa": len(sequence),
            "cxxch_count": len(motifs),
            "cxxch_positions": ";".join(str(position) for position, _ in motifs),
            "cxxch_motifs": ";".join(motif for _, motif in motifs),
        })

    carrier_mags = {(row["dataset_id"], row["genome_id"]) for row in rows}
    carrier_ani = {row["ani_cluster_95"] for row in rows}
    if len(carrier_mags) != 19 or len(carrier_ani) != 19:
        raise RuntimeError(
            f"Frozen carrier gate failed: MAGs={len(carrier_mags)}, ANI95={len(carrier_ani)}"
        )

    args.outdir.mkdir(parents=True)
    fasta = args.outdir / "DPFQ008_family_20.faa"
    with fasta.open("x", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f">{row['stable_id']} {row['protein_id']}\n")
            sequence = sequences[row["protein_id"]]
            handle.write("\n".join(sequence[i:i + 80] for i in range(0, len(sequence), 80)) + "\n")

    table = args.outdir / "DPFQ008_family_architecture.tsv"
    with table.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    summary = [
        "DPFQ008 frozen-family architecture audit",
        f"family_members={len(rows)}",
        f"carrier_MAGs={len(carrier_mags)}",
        f"carrier_ANI95_clusters={len(carrier_ani)}",
        f"members_with_2plus_CXXCH={sum(int(row['cxxch_count']) >= 2 for row in rows)}",
        f"members_with_any_CXXCH={sum(int(row['cxxch_count']) >= 1 for row in rows)}",
        f"length_min={min(int(row['length_aa']) for row in rows)}",
        f"length_median={sorted(int(row['length_aa']) for row in rows)[len(rows)//2]}",
        f"length_max={max(int(row['length_aa']) for row in rows)}",
        f"family_fasta_sha256={digest(fasta)}",
        "claim_boundary=motif conservation supports haem-binding architecture but does not establish substrate, flux, or activity",
        "status=PASS",
    ]
    (args.outdir / "run_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("\n".join(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
