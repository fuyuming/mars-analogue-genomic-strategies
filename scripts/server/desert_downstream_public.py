#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Public-facing, path-parameterised form of the round-39 dryland downstream
pipeline (CheckM2 quality gate -> primary cohort -> Prodigal/KOfam annotation
-> GTDB-Tk taxonomy). No ecological model is fitted.

INTEGRATION-NOT-RETESTED: faithful de-privatised transcription of the authoring
script. It has NOT been run end-to-end here (needs CheckM2 v1.1.0, the KOfamScan
1.3.0 profiles, GTDB-Tk r226 and the Prodigal outputs on a compute server).
Archived analysis code, not a verified smoke test.

All tool/database locations are command-line arguments; no private host path is
embedded. The project directory must contain `prediction_summary.json`
(verified/failed counts), `manifest.tsv` (columns: accession) and
`prodigal/<accession>/{genome.fna,prefixed_min20.faa}`.

Example
-------
python desert_downstream_public.py \
    --workdir   /path/to/desert_project \
    --checkm2-env /path/to/envs/checkm2-1.1.0 \
    --checkm2-db  /path/to/CheckM2_database/uniref100.KO.1.dmnd \
    --kofam-dir   /path/to/kofam_scan-1.3.0 \
    --gtdb-bin    /path/to/envs/gtdbtk-2.5.2/bin \
    --gtdb-db     /path/to/gtdb_r226 \
    --threads     48
"""
import argparse
import concurrent.futures
import csv
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(description="Dryland quality->cohort->annotation pipeline (archived, not retested)")
    ap.add_argument("--workdir", type=Path, required=True, help="project directory (manifest.tsv, prodigal/)")
    ap.add_argument("--checkm2-env", type=Path, required=True, help="env dir containing bin/checkm2")
    ap.add_argument("--checkm2-db", type=Path, required=True, help="CheckM2 diamond database")
    ap.add_argument("--kofam-dir", type=Path, required=True, help="KofamScan 1.3.0 directory (exec_annotation, config.yml)")
    ap.add_argument("--gtdb-bin", type=Path, required=True, help="directory containing the gtdbtk executable")
    ap.add_argument("--gtdb-db", type=Path, required=True, help="GTDB-Tk r226 full_package reference")
    ap.add_argument("--threads", type=int, default=48)
    ap.add_argument("--expect-genomes", type=int, default=911)
    return ap.parse_args()


def main():
    a = parse_args()
    R = a.workdir.resolve()
    CE = a.checkm2_env.resolve()
    DB = a.checkm2_db.resolve()
    K = a.kofam_dir.resolve()
    GE = a.gtdb_bin.resolve()
    GDB = a.gtdb_db.resolve()

    def run(cmd, log, env=None):
        with log.open("a") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S") + " CMD " + repr(cmd) + "\n")
            f.flush()
            subprocess.run(list(map(str, cmd)), stdout=f, stderr=f, check=True, env=env)

    def status(stage, **kwargs):
        (R / "downstream_status.json").write_text(json.dumps(
            {"stage": stage, "time": time.strftime("%Y-%m-%d %H:%M:%S"), **kwargs}, indent=2))

    s = json.loads((R / "prediction_summary.json").read_text())
    assert s["verified"] == a.expect_genomes and s["failed"] == 0
    manifest = list(csv.DictReader((R / "manifest.tsv").open(), delimiter="\t"))
    inputs = R / "quality_inputs"
    inputs.mkdir(exist_ok=True)
    for r in manifest:
        src = R / "prodigal" / r["accession"] / "genome.fna"
        dst = inputs / (r["accession"] + ".fna")
        assert src.exists()
        if not dst.exists():
            dst.symlink_to(src)

    env = os.environ.copy()
    env["PATH"] = str(CE / "bin") + os.pathsep + env["PATH"]
    qc = R / "checkm2_1.1.0"
    report = qc / "quality_report.tsv"
    status("CheckM2_running", n_genomes=a.expect_genomes, threads=a.threads)
    if not report.exists():
        run([CE / "bin/checkm2", "predict", "--input", inputs, "--output-directory", qc,
             "--threads", a.threads, "--extension", "fna", "--database_path", DB],
            R / "checkm2_wrapper.log", env)
    rows = list(csv.DictReader(report.open(), delimiter="\t"))
    assert len(rows) == a.expect_genomes and len({r["Name"] for r in rows}) == a.expect_genomes
    accepted = [r for r in rows if float(r["Completeness"]) >= 50 and float(r["Contamination"]) < 10]
    strict = [r for r in accepted if float(r["Completeness"]) >= 90 and float(r["Contamination"]) < 5]
    (R / "quality_gate_summary.json").write_text(json.dumps(
        {"total": a.expect_genomes, "primary_50_10": len(accepted), "strict_90_5": len(strict),
         "excluded": a.expect_genomes - len(accepted),
         "thresholds": "completeness>=50,contamination<10; strict>=90,<5; not automatically paperMIMAGHQ"}, indent=2))
    with (R / "cohort_quality_registry.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]) + ["primary", "strict"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({**r, "primary": r in accepted, "strict": r in strict})

    primary = R / "primary_inputs"
    primary.mkdir(exist_ok=True)
    combined = R / "desert_primary_uniform.proteins.faa"
    nproteins = 0
    with combined.with_suffix(".part").open("wb") as out:
        for r in sorted(accepted, key=lambda r: r["Name"]):
            name = r["Name"]
            src = inputs / (name + ".fna")
            dst = primary / (name + ".fna")
            if not dst.exists():
                dst.symlink_to(src)
            data = (R / "prodigal" / name / "prefixed_min20.faa").read_bytes()
            nproteins += sum(l.startswith(b">") for l in data.splitlines())
            out.write(data)
    combined.with_suffix(".part").rename(combined)
    status("annotation_and_taxonomy_running", primary_MAGs=len(accepted), strict_MAGs=len(strict), proteins=nproteins)

    meta = {
        "quality_report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
        "combined_proteins_sha256": hashlib.sha256(combined.read_bytes()).hexdigest(),
        "kofam_config_sha256": hashlib.sha256((K / "config.yml").read_bytes()).hexdigest(),
        "ko_list_sha256": hashlib.sha256((K.parent / "ko_list").read_bytes()).hexdigest(),
        "prokaryote_hal_sha256": hashlib.sha256((K.parent / "profiles/prokaryote.hal").read_bytes()).hexdigest(),
        "gtdb_metadata_sha256": hashlib.sha256((GDB / "metadata/metadata.txt").read_bytes()).hexdigest(),
        "checkm2_database_bytes": DB.stat().st_size, "full_profiles": True,
    }
    (R / "annotation_provenance.json").write_text(json.dumps(meta, indent=2))

    def kofam():
        out = R / "full_kofam_1.3.0"
        out.mkdir(exist_ok=True)
        final = out / "mapper_one_line.tsv"
        done = out / "verified.json"
        if done.exists():
            return
        assert not final.exists(), "Unverified existing result: inspect before restart"
        run([K / "exec_annotation", "-c", K / "config.yml", "--cpu", a.threads, "--tmp-dir", out / "tmp",
             "--format", "mapper-one-line", "--no-report-unannotated", "-o", final, combined], out / "run.log")
        seen, bad = set(), 0
        with final.open() as f:
            for l in f:
                if not l.strip() or l.startswith("#"):
                    continue
                parts = l.rstrip().split("\t")
                seen.add(parts[0])
                bad += not parts[0].startswith("GCA_")
        assert bad == 0 and len(seen) <= nproteins
        done.write_text(json.dumps({"annotated_proteins": len(seen), "input_proteins": nproteins,
                                    "malformed": bad, "sha256": hashlib.sha256(final.read_bytes()).hexdigest()}, indent=2))

    def gtdb():
        out = R / "gtdbtk_r226"
        env2 = os.environ.copy()
        env2["PATH"] = str(GE) + os.pathsep + env2["PATH"]
        env2["GTDBTK_DATA_PATH"] = str(GDB)
        done = R / "gtdb_verified.json"
        if done.exists():
            return
        if out.exists():
            raise RuntimeError("Inspect incomplete GTDB output; no implicit overwrite")
        run([GE / "gtdbtk", "classify_wf", "--genome_dir", primary, "--out_dir", out,
             "--extension", "fna", "--cpus", a.threads, "--pplacer_cpus", str(a.threads // 2)],
            R / "gtdb_wrapper.log", env2)
        ids = []
        for name in ["gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"]:
            p = out / "classify" / name
            if p.exists():
                ids += [r["user_genome"] for r in csv.DictReader(p.open(), delimiter="\t")]
        fail = out / "identify/gtdbtk.failed_genomes.tsv"
        if fail.exists():
            ids += [r["genome_id"] for r in csv.DictReader(fail.open(), delimiter="\t")]
        assert len(ids) == len(accepted) and len(set(ids)) == len(ids)
        done.write_text(json.dumps({"accounted": len(ids), "expected": len(accepted)}, indent=2))

    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        for f in [pool.submit(kofam), pool.submit(gtdb)]:
            f.result()
    status("annotation_and_taxonomy_complete", primary_MAGs=len(accepted),
           note="ANI deduplication and ecological modelling still pending; no new discovery claimed")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:  # noqa: BLE001
        # best-effort status write using the default workdir supplied on the CLI
        try:
            import sys
            wd = Path(sys.argv[sys.argv.index("--workdir") + 1])
            (wd / "downstream_status.json").write_text(json.dumps({"stage": "failed", "error": repr(e)}, indent=2))
        except Exception:  # noqa: BLE001
            pass
        raise
