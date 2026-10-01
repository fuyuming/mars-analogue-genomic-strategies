#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Public-facing, path-parameterised form of the round-53 strict-cohort
marker-extraction + tree-inference driver (GTDB-Tk r226 markers -> IQ-TREE).

INTEGRATION-NOT-RETESTED: this wrapper is a faithful, de-privatised transcription
of the authoring script. It has NOT been executed end-to-end here, because it
requires a GTDB-Tk environment, the r226 reference package, IQ-TREE2 and the
1129-genome strict manifest on a compute server. Treat it as archived analysis
code, not as a verified smoke test.

Every filesystem location is supplied on the command line; no private host path
is embedded.

Example
-------
python run_joint_strict_public.py \
    --workdir  /path/to/joint_strict_project \
    --gtdb-bin /path/to/envs/gtdbtk-2.5.2/bin \
    --gtdb-db  /path/to/gtdb_r226 \
    --iqtree2  /usr/bin/iqtree2 \
    --threads  32 \
    --expect-tips 1129

The project directory must already contain `manifest.tsv` (columns: tip_id,
genome_id, domain), `preflight.json` (key: input_manifest_sha256) and
`genomes/<accession>.fna`.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
import re
import subprocess
import time
import traceback
from pathlib import Path


def parse_args():
    ap = argparse.ArgumentParser(description="Strict-cohort GTDB-Tk marker + IQ-TREE2 pipeline (archived, not retested)")
    ap.add_argument("--workdir", type=Path, required=True, help="project directory (manifest.tsv, genomes/)")
    ap.add_argument("--gtdb-bin", type=Path, required=True, help="directory containing the gtdbtk executable")
    ap.add_argument("--gtdb-db", type=Path, required=True, help="GTDB-Tk reference package (release226 full_package)")
    ap.add_argument("--iqtree2", type=Path, default=Path("iqtree2"), help="IQ-TREE2 executable (default: iqtree2 on PATH)")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--expect-tips", type=int, default=1129, help="expected number of strict tips")
    ap.add_argument("--seed", type=int, default=20261001)
    return ap.parse_args()


def main():
    a = parse_args()
    R = a.workdir.resolve()
    E = a.gtdb_bin.resolve()
    DB = a.gtdb_db.resolve()

    lock = (R / "pipeline.lock").open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    rows = list(csv.DictReader((R / "manifest.tsv").open(), delimiter="\t"))
    expected = {x["tip_id"]: x for x in rows}
    manifest_hash = hashlib.sha256((R / "manifest.tsv").read_bytes()).hexdigest()
    assert len(expected) == a.expect_tips, (len(expected), a.expect_tips)
    assert json.loads((R / "preflight.json").read_text())["input_manifest_sha256"] == manifest_hash

    env = os.environ.copy()
    env["PATH"] = str(E) + os.pathsep + env["PATH"]
    env["GTDBTK_DATA_PATH"] = str(DB)

    def state(stage, **kwargs):
        d = {"stage": stage, "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
             "pid": os.getpid(), "manifest_sha256": manifest_hash, **kwargs}
        p = R / "status.json.tmp"
        p.write_text(json.dumps(d, indent=2))
        p.replace(R / "status.json")
        print(json.dumps(d), flush=True)

    def run(name, args):
        receipt = R / (name + ".done.json")
        if receipt.exists():
            done = json.loads(receipt.read_text())
            assert done["args"] == list(map(str, args)) and done["manifest_sha256"] == manifest_hash
            return
        state(name, args=list(map(str, args)), threads=a.threads)
        with (R / (name + ".log")).open("a") as f:
            t = time.time()
            p = subprocess.run(list(map(str, args)), cwd=R, env=env, stdout=f, stderr=subprocess.STDOUT)
        assert p.returncode == 0, (name, p.returncode)
        receipt.write_text(json.dumps({"args": list(map(str, args)), "manifest_sha256": manifest_hash,
                                       "seconds": time.time() - t, "exit": p.returncode}, indent=2))

    def fasta(path):
        seqs, key = {}, None
        with gzip.open(path, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    key = line[1:].split()[0]
                    assert key not in seqs
                    seqs[key] = ""
                else:
                    seqs[key] += line.strip()
        return seqs

    try:
        run("identify", [E / "gtdbtk", "identify", "--genome_dir", R / "genomes",
                         "--out_dir", R / "markers", "--extension", "fna", "--cpus", a.threads])
        run("align", [E / "gtdbtk", "align", "--identify_dir", R / "markers",
                      "--out_dir", R / "alignment", "--skip_gtdb_refs", "--cpus", a.threads])
        summaries, allids = [], set()
        out = R / "inference"
        out.mkdir(exist_ok=True)
        for marker, domain, n, L in [("bac120", "Bacteria", 1092, 5036), ("ar53", "Archaea", 37, 8062)]:
            seqs = fasta(R / f"alignment/align/gtdbtk.{marker}.user_msa.fasta.gz")
            want = {k for k, x in expected.items() if x["domain"] == domain}
            missing, extra = sorted(want - set(seqs)), sorted(set(seqs) - want)
            summaries.append({"marker": marker, "expected": n, "observed": len(seqs),
                              "missing": missing, "unexpected": extra,
                              "lengths": sorted(set(map(len, seqs.values())))})
            (R / "alignment_qc.json").write_text(json.dumps(summaries, indent=2))
            assert not missing and not extra and len(seqs) == n, (marker, "tip/domain mismatch", missing, extra)
            assert set(map(len, seqs.values())) == {L}
            assert not allids & set(seqs)
            allids |= set(seqs)
            records = []
            for k, s in seqs.items():
                occupancy = sum(c not in "-.X?" for c in s) / L
                assert occupancy >= .10
                records.append({"tip_id": k, "genome_id": expected[k]["genome_id"], "length": L,
                                "observed_aa_fraction": occupancy})
            with (R / f"{marker}_occupancy.tsv").open("w", newline="") as f:
                w = csv.DictWriter(f, records[0].keys(), delimiter="\t")
                w.writeheader()
                w.writerows(records)
            content = "".join(">" + k + "\n" + seqs[k] + "\n" for k in sorted(seqs))
            p = out / f"{marker}_strict.fasta"
            if p.exists():
                assert p.read_text() == content
            else:
                p.write_text(content)
        assert allids == set(expected)
        state("alignment_qc_passed", domains=summaries)
        for marker in ["bac120", "ar53"]:
            prefix = out / f"{marker}_strict_LG_F_R10"
            run("iqtree_" + marker, [a.iqtree2, "-s", out / f"{marker}_strict.fasta", "-m", "LG+F+R10",
                                     "-B", "1000", "--alrt", "1000", "-T", a.threads, "--mem", "64G",
                                     "-seed", a.seed, "--prefix", prefix])
            tree = Path(str(prefix) + ".treefile")
            assert tree.exists()
            tips = re.findall(r"\bNEE53_S_\d{5}\b", tree.read_text())
            seqids = set(fasta(R / f"alignment/align/gtdbtk.{marker}.user_msa.fasta.gz"))
            assert len(tips) == len(seqids) and set(tips) == seqids
        state("complete", trees=[str(out / f"{m}_strict_LG_F_R10.treefile") for m in ["bac120", "ar53"]],
              note="Tree inference complete; ecological comparisons and topology sensitivity not yet run")
    except Exception as ex:  # noqa: BLE001
        state("failed", error=repr(ex))
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
