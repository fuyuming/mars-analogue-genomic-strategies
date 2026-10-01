#!/usr/bin/env python3
"""Matched within-MAG permutation test for DPFQ008 cytochrome/haem context.

Controls are proteins from the supplied sequence universe and the same carrier
MAG, matched to each target on length (+/-20%) and complete versus truncated
+/-5-gene window. Controls overlapping a DPFQ008 target window are excluded.
The primary statistic is the number of carrier MAGs with at least one positive
target context; permutations sample one matched control per target stratum and
aggregate duplicated targets at MAG level. Frozen MMseqs-cluster and expanded
HMM-hit target sets are both supported through explicit gates.
"""

from __future__ import annotations

import argparse
import csv
import random
import re
from collections import Counter, defaultdict
from pathlib import Path


REPRESENTATIVE = (
    "qaidam_basin_mmag_1773|qaidam_basin_mmag_1773_CHGRbin21|"
    "k141_121601_length_11112_cov_23.0792_5"
)
EXPECTED_TARGETS = 20
CONTEXT_RE = re.compile(r"cytochrome|CXXCH|heme|haem|porphyr", re.I)


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--cluster-tsv", type=Path)
    source.add_argument("--hmm-hit-table", type=Path)
    proteins = p.add_mutually_exclusive_group(required=True)
    proteins.add_argument("--unknown-fasta", type=Path)
    proteins.add_argument("--protein-fasta", type=Path)
    p.add_argument("--prodigal-root", required=True, type=Path)
    p.add_argument(
        "--prodigal-dataset-root", action="append", default=[],
        help="Optional DATASET=/path/to/dataset-directory override; repeatable",
    )
    p.add_argument("--eggnog-root", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--flank", type=int, default=5)
    p.add_argument("--length-tolerance", type=float, default=0.20)
    p.add_argument("--permutations", type=int, default=100000)
    p.add_argument("--seed", type=int, default=20260719)
    p.add_argument("--minimum-controls", type=int, default=10)
    p.add_argument("--hmm-threshold", default="strict")
    p.add_argument("--expected-targets", type=int, default=EXPECTED_TARGETS)
    p.add_argument("--expected-carrier-mags", type=int, default=19)
    return p.parse_args()


def split_id(value: str) -> tuple[str, str, str]:
    fields = value.split("|", 2)
    if len(fields) != 3:
        raise ValueError(value)
    return fields[0], fields[1], fields[2]


def targets(path: Path, expected: int) -> list[str]:
    result = []
    with path.open(encoding="utf-8") as handle:
        for i, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2:
                raise ValueError(f"Malformed cluster row {i}")
            if fields[0] == REPRESENTATIVE:
                result.append(fields[1])
    if len(result) != expected or len(set(result)) != expected:
        raise RuntimeError(f"Frozen target gate failed: {len(result)}")
    return result


def hmm_targets(path: Path, threshold: str, expected: int) -> list[str]:
    with path.open(encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if threshold not in (reader.fieldnames or []):
            raise RuntimeError(f"HMM threshold column absent: {threshold}")
        result = [row["protein_id"] for row in reader if row[threshold] == "1"]
    if len(result) != expected or len(set(result)) != expected:
        raise RuntimeError(f"Frozen HMM target gate failed: {len(result)}")
    return result


def unresolved_in_carriers(path: Path, carrier_keys: set[tuple[str, str]]) -> dict[tuple[str, str], dict[str, int]]:
    result: dict[tuple[str, str], dict[str, int]] = defaultdict(dict)
    current: tuple[tuple[str, str], str] | None = None
    length = 0

    def commit() -> None:
        if current is not None and length > 0:
            key, original = current
            if original in result[key] and result[key][original] != length:
                raise RuntimeError(f"Conflicting FASTA duplicate: {key}|{original}")
            result[key][original] = length

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                commit()
                full = line[1:].split()[0]
                dataset, genome, original = split_id(full)
                key = (dataset, genome)
                current = (key, original) if key in carrier_keys else None
                length = 0
            elif current is not None:
                length += len(line.strip().rstrip("*"))
        commit()
    return result


def faa_id_map(path: Path) -> dict[str, str]:
    mapping = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            original = line[1:].split()[0]
            match = re.search(r"\bID=([^;\s]+)", line)
            if match:
                mapping[original] = match.group(1)
    return mapping


def gff_rows(path: Path) -> tuple[list[dict[str, str]], dict[str, tuple[str, int]], dict[str, dict[str, str]]]:
    by_contig: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) != 9 or f[2] != "CDS":
                continue
            attrs = dict(token.split("=", 1) for token in f[8].split(";") if "=" in token)
            by_contig[f[0]].append({
                "seqid": f[0], "gff_id": attrs.get("ID", ""),
                "start": f[3], "end": f[4], "strand": f[6],
            })
    flat = []
    locations = {}
    for seqid, rows in by_contig.items():
        for index, row in enumerate(rows):
            locations[row["gff_id"]] = (seqid, index)
            flat.append(row)
    details = {row["gff_id"]: row for row in flat}
    return flat, locations, details


def original_query(dataset: str, genome: str, seqid: str, gff_id: str) -> str:
    local = gff_id.rsplit("_", 1)[1]
    return f"{dataset}|{genome}|{seqid}_{local}"


def eggnog_path(root: Path, dataset: str) -> Path:
    matches = sorted(root.glob(f"{dataset}*_dbmem.emapper.annotations"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one eggNOG file for {dataset}, saw {matches}")
    return matches[0]


def main() -> int:
    a = args()
    target_path = a.cluster_tsv or a.hmm_hit_table
    sequence_fasta = a.unknown_fasta or a.protein_fasta
    for path in (target_path, sequence_fasta, a.prodigal_root, a.eggnog_root):
        path.resolve(strict=True)
    if a.outdir.exists():
        raise FileExistsError(a.outdir)
    if not (0 < a.length_tolerance < 1) or a.flank < 1 or a.permutations < 999:
        raise ValueError("Invalid matching/permutation parameters")
    prodigal_overrides = {}
    for item in a.prodigal_dataset_root:
        if "=" not in item:
            raise ValueError(f"Malformed Prodigal dataset override: {item}")
        dataset, value = item.split("=", 1)
        path = Path(value)
        path.resolve(strict=True)
        if dataset in prodigal_overrides:
            raise ValueError(f"Duplicate Prodigal dataset override: {dataset}")
        prodigal_overrides[dataset] = path

    target_ids = (
        targets(a.cluster_tsv, a.expected_targets) if a.cluster_tsv else
        hmm_targets(a.hmm_hit_table, a.hmm_threshold, a.expected_targets)
    )
    target_set = set(target_ids)
    carrier_keys = {split_id(value)[:2] for value in target_ids}
    if len(carrier_keys) != a.expected_carrier_mags:
        raise RuntimeError(f"Expected {a.expected_carrier_mags} carrier MAGs, saw {len(carrier_keys)}")
    unresolved = unresolved_in_carriers(sequence_fasta, carrier_keys)

    records = {}
    neighborhoods: dict[str, list[str]] = {}
    target_windows: dict[tuple[str, str], set[tuple[str, int]]] = defaultdict(set)
    for dataset, genome in sorted(carrier_keys):
        dataset_root = prodigal_overrides.get(dataset, a.prodigal_root / dataset)
        faa = dataset_root / f"{genome}.faa"
        gff = dataset_root / f"{genome}.gff"
        mapping = faa_id_map(faa)
        _, locations, details = gff_rows(gff)
        contig_rows: dict[str, list[str]] = defaultdict(list)
        for gff_id, (seqid, index) in locations.items():
            if len(contig_rows[seqid]) <= index:
                contig_rows[seqid].extend([""] * (index + 1 - len(contig_rows[seqid])))
            contig_rows[seqid][index] = gff_id
        for original, length in unresolved[(dataset, genome)].items():
            gff_id = mapping.get(original)
            if not gff_id or gff_id not in locations:
                continue
            seqid, index = locations[gff_id]
            genes = contig_rows[seqid]
            left, right = min(index, a.flank), min(len(genes) - index - 1, a.flank)
            full = left == a.flank and right == a.flank
            full_id = f"{dataset}|{genome}|{original}"
            neighbor_ids = [
                original_query(dataset, genome, seqid, genes[j])
                for j in range(max(0, index-a.flank), min(len(genes), index+a.flank+1))
                if j != index and genes[j]
            ]
            records[full_id] = {
                "dataset": dataset, "genome": genome, "original": original,
                "length": length, "seqid": seqid, "index": index,
                "full_window": full, "strand": details[gff_id]["strand"],
                "start": int(details[gff_id]["start"]), "end": int(details[gff_id]["end"]),
            }
            neighborhoods[full_id] = neighbor_ids
        for target in [x for x in target_ids if split_id(x)[:2] == (dataset, genome)]:
            r = records.get(target)
            if r is None:
                raise RuntimeError(f"Target not resolved in Prodigal files: {target}")
            for j in range(max(0, r["index"]-a.flank), r["index"]+a.flank+1):
                target_windows[(dataset, genome)].add((r["seqid"], j))

    pools = {}
    for target in target_ids:
        tr = records[target]
        low = tr["length"] * (1-a.length_tolerance)
        high = tr["length"] * (1+a.length_tolerance)
        key = (tr["dataset"], tr["genome"])
        eligible = []
        for protein, row in records.items():
            if protein in target_set or (row["dataset"], row["genome"]) != key:
                continue
            if not low <= row["length"] <= high or row["full_window"] != tr["full_window"]:
                continue
            if (row["seqid"], row["index"]) in target_windows[key]:
                continue
            eligible.append(protein)
        if len(eligible) < a.minimum_controls:
            raise RuntimeError(f"Too few matched controls for {target}: {len(eligible)}")
        pools[target] = sorted(eligible)

    required = {query for protein in set(target_ids).union(*map(set, pools.values())) for query in neighborhoods[protein]}
    annotation_text = {}
    annotations = {}
    for dataset in {key[0] for key in carrier_keys}:
        with eggnog_path(a.eggnog_root, dataset).open(errors="replace", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[0] in required:
                    annotation_text[fields[0]] = " ".join((fields[7], fields[8], fields[20]))
                    annotations[fields[0]] = {
                        "Description": fields[7], "Preferred_name": fields[8],
                        "KEGG_ko": fields[11], "PFAMs": fields[20],
                    }

    positive = {
        protein: any(CONTEXT_RE.search(annotation_text.get(query, "")) for query in neighborhoods[protein])
        for protein in set(target_ids).union(*map(set, pools.values()))
    }
    targets_by_mag: dict[tuple[str, str], list[str]] = defaultdict(list)
    for target in target_ids:
        tr = records[target]
        targets_by_mag[(tr["dataset"], tr["genome"])].append(target)
    observed_mag = sum(any(positive[t] for t in ts) for ts in targets_by_mag.values())
    observed_protein = sum(positive[t] for t in target_ids)

    rng = random.Random(a.seed)
    null_mag = []
    null_protein = []
    for _ in range(a.permutations):
        sampled = {target: rng.choice(pools[target]) for target in target_ids}
        null_protein.append(sum(positive[value] for value in sampled.values()))
        null_mag.append(sum(any(positive[sampled[t]] for t in ts) for ts in targets_by_mag.values()))
    p_mag = (1 + sum(value >= observed_mag for value in null_mag)) / (a.permutations + 1)
    p_protein = (1 + sum(value >= observed_protein for value in null_protein)) / (a.permutations + 1)

    def quantile(values: list[int], p: float) -> int:
        ordered = sorted(values)
        return ordered[int((len(ordered)-1)*p)]

    a.outdir.mkdir(parents=True)
    with (a.outdir / "matched_target_control_strata.tsv").open("x", newline="", encoding="utf-8") as handle:
        fields = ["target_protein_id", "genome_id", "target_length", "target_full_window",
                  "target_context_positive", "matched_controls", "positive_controls", "positive_control_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for target in sorted(target_ids):
            tr = records[target]
            hits = sum(positive[x] for x in pools[target])
            writer.writerow({
                "target_protein_id": target, "genome_id": tr["genome"],
                "target_length": tr["length"], "target_full_window": int(tr["full_window"]),
                "target_context_positive": int(positive[target]), "matched_controls": len(pools[target]),
                "positive_controls": hits, "positive_control_fraction": f"{hits/len(pools[target]):.6f}",
            })
    summary_rows = [
        ("MAG", observed_mag, len(targets_by_mag), sum(null_mag)/len(null_mag), quantile(null_mag,.025), quantile(null_mag,.975), p_mag),
        ("protein", observed_protein, len(target_ids), sum(null_protein)/len(null_protein), quantile(null_protein,.025), quantile(null_protein,.975), p_protein),
    ]
    with (a.outdir / "matched_permutation_summary.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["level","observed_positive","denominator","null_mean","null_q025","null_q975","empirical_p_ge"])
        writer.writerows(summary_rows)
    with (a.outdir / "null_distribution_frequency.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["level", "positive_count", "permutation_frequency", "permutations"])
        for level, values in (("MAG", null_mag), ("protein", null_protein)):
            for count, frequency in sorted(Counter(values).items()):
                writer.writerow([level, count, frequency, len(values)])

    neighborhood_fields = [
        "target_protein_id", "dataset_id", "genome_id", "target_full_window",
        "target_context_positive", "neighbor_protein_id", "relative_gene",
        "oriented_gene", "neighbor_strand", "oriented_strand", "gene_class",
        "Description", "Preferred_name", "KEGG_ko", "PFAMs",
    ]
    with (a.outdir / "expanded_target_neighborhoods.tsv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=neighborhood_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for target in sorted(target_ids):
            tr = records[target]
            orient = 1 if tr["strand"] == "+" else -1
            members = [(target, 0)]
            members.extend((query, records[query]["index"] - tr["index"]) for query in neighborhoods[target])
            for query, relative in sorted(members, key=lambda item: item[1]):
                row = records[query]
                ann = annotations.get(query, {"Description": "-", "Preferred_name": "-", "KEGG_ko": "-", "PFAMs": "-"})
                text = " ".join(ann.values())
                gene_class = (
                    "target" if query == target else
                    "cyto_b_N" if "Cytochrom_B_N_2" in ann["PFAMs"] else
                    "cyto_b_C" if "Cytochrom_B_C" in ann["PFAMs"] else
                    "heme" if CONTEXT_RE.search(text) else
                    "other" if any(value not in ("", "-") for value in ann.values()) else
                    "unresolved"
                )
                oriented_strand = row["strand"] if orient == 1 else ("-" if row["strand"] == "+" else "+")
                writer.writerow({
                    "target_protein_id": target, "dataset_id": tr["dataset"], "genome_id": tr["genome"],
                    "target_full_window": int(tr["full_window"]),
                    "target_context_positive": int(positive[target]), "neighbor_protein_id": query,
                    "relative_gene": relative, "oriented_gene": relative * orient,
                    "neighbor_strand": row["strand"], "oriented_strand": oriented_strand,
                    "gene_class": gene_class, **ann,
                })
    lines = [
        "DPFQ008 matched within-MAG context enrichment",
        f"targets={len(target_ids)}", f"carrier_MAGs={len(targets_by_mag)}",
        f"length_tolerance={a.length_tolerance}", f"flank={a.flank}",
        f"permutations={a.permutations}", f"seed={a.seed}",
        f"observed_positive_MAGs={observed_mag}/{len(targets_by_mag)}", f"empirical_MAG_p_ge={p_mag:.8g}",
        f"observed_positive_proteins={observed_protein}/{len(target_ids)}", f"empirical_protein_p_ge={p_protein:.8g}",
        "control_scope=same MAG; input protein universe; matched length and complete/truncated window; target-window overlaps excluded",
        "claim_boundary=enrichment supports non-random genomic context but not substrate, flux, expression, or activity",
        "status=PASS",
    ]
    (a.outdir / "run_summary.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
