"""Portable CLI: descriptive joint-configuration support diagnostic.

Publication-facing, path-parameterised rendering of the frozen analysis in
``analysis/analysis_contract.md`` (see ``data/upstream/joint/analysis_contract.md``).
No absolute private path is present anywhere in this file: every input location
is resolved from ``--data-root`` (default: this repository's ``data/upstream``).

Frozen scientific rules (DO NOT alter):
  * four pre-specified maintenance<->Cox/Cut pairs (KdpABC, ProVWX, EctABC, BetAB
    each with K03518/K03519/K03520);
  * Cox/Cut co-detection is NOT a confirmed CO-oxidation marker;
  * joint-state support gate = every one of the four binary states has >= min
    representative MAGs within a family AND >= 2 source/site groups each with
    >= 3 MAGs; thresholds 2 / 3 / 5 are reported as declared sensitivities;
  * no p-values, no environment model, no causal reading.

Outputs (written to --outdir): joint_summary.tsv, family_joint_cells.tsv,
joint_support.tsv, genome_system_states.tsv, validation.json.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

RULES = {
    "KdpABC": ["K01546", "K01547", "K01548"],
    "ProVWX": ["K02000", "K02001", "K02002"],
    "EctABC": ["K06718", "K00836", "K06720"],
    "BetAB": ["K00108", "K00130"],
    "CoxCut": ["K03518", "K03519", "K03520"],
}


def run(data_root: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    states_dir = data_root / "cohorts"
    original_presence = data_root / "panel_repair/amended_inputs/restored74_presence.csv.gz"
    new_presence = data_root / "cohorts/new_catalogue/restored74_presence.tsv.gz"

    raw = pd.read_csv(original_presence)
    tax = raw[["genome_id", "family", "genus"]].drop_duplicates().set_index("genome_id")
    assert not tax.index.duplicated().any()
    old = raw.pivot(index="genome_id", columns="KO", values="present")
    new = pd.read_csv(new_presence, sep="\t", index_col=0)

    summaries, cells, support, allstates = [], [], [], []
    for cohort in ["original_primary", "original_strict", "new_desert_primary", "new_desert_strict"]:
        states = pd.read_csv(states_dir / f"{cohort}_system_states.tsv", sep="\t").set_index("genome_id")
        if cohort.startswith("original"):
            meta = states.join(tax)
            meta["unit"] = meta.dataset_id
            mat = old.loc[meta.index]
        else:
            suffix = cohort.replace("new_desert_", "")
            meta = pd.read_csv(states_dir / f"new_desert_{suffix}_independent_representatives.tsv",
                               sep="\t").set_index("genome_id")
            meta = meta.loc[states.index]
            meta["unit"] = meta.site_labels
            mat = new.loc[meta.index]
        assert mat.index.is_unique and mat.loc[:, sum(RULES.values(), [])].notna().all().all()

        fam = meta.family.fillna("")
        known = ~fam.isin(["", "f__", "unclassified", "unassigned", "nan"])
        fam = fam.where(known)

        calls = {}
        for name, kos in RULES.items():
            hits = mat[kos].sum(axis=1)
            calls[name] = hits.eq(len(kos))
            c = hits.map(lambda n: "all_detected" if n == len(kos) else ("partial" if n > 0 else "not_detected"))
            if name != "CoxCut":
                assert (calls[name] == states[name].eq("panel_complete")).all(), "prior maintenance calls changed"
            for gid in meta.index:
                allstates.append({"cohort": cohort, "genome_id": gid, "family": fam.loc[gid],
                                  "unit": meta.loc[gid, "unit"], "system": name, "state": c.loc[gid]})

        for name in list(RULES)[:-1]:
            d = pd.DataFrame({"maintenance": calls[name].astype(int), "energy": calls["CoxCut"].astype(int),
                              "family": fam, "unit": meta.unit})
            d["joint"] = 2 * d.maintenance + d.energy
            ct = d.joint.value_counts()
            eligible, familyexp, familyobs, nk = [], 0, 0, 0
            for family, g in d.dropna(subset=["family"]).groupby("family"):
                counts = [int((g.joint == j).sum()) for j in range(4)]
                units = g.unit.value_counts()
                n = len(g)
                exp = float(g.maintenance.sum() * g.energy.sum() / n)
                familyexp += exp
                familyobs += counts[3]
                nk += n
                cells.append({"cohort": cohort, "pair": name + "__CoxCut", "family": family, "n": n,
                              "n00": counts[0], "n01": counts[1], "n10": counts[2], "n11": counts[3],
                              "source_or_site_groups": len(units), "groups_ge3": int((units >= 3).sum()),
                              "expected_both_within_family": exp})
                for gate in [2, 3, 5]:
                    if min(counts) >= gate and (units >= 3).sum() >= 2:
                        eligible.append((gate, family, n, len(units)))
            for gate in [2, 3, 5]:
                ee = [x for x in eligible if x[0] == gate]
                support.append({"cohort": cohort, "pair": name + "__CoxCut", "min_per_joint_state": gate,
                                "eligible_families": len(ee), "eligible_genomes": sum(x[2] for x in ee),
                                "families": ";".join(x[1] for x in ee)})
            fam_d = d.dropna(subset=["family"])
            exp_pooled_assigned = float(fam_d.maintenance.sum() * fam_d.energy.sum() / len(fam_d))
            summaries.append({"cohort": cohort, "pair": name + "__CoxCut", "n": len(d),
                              "n00": int(ct.get(0, 0)), "n01": int(ct.get(1, 0)),
                              "n10": int(ct.get(2, 0)), "n11": int(ct.get(3, 0)),
                              "expected_both_pooled": float(d.maintenance.sum() * d.energy.sum() / len(d)),
                              "family_assigned_n": nk, "both_family_assigned": familyobs,
                              "expected_both_family_conditioned": familyexp,
                              "expected_both_pooled_family_assigned": exp_pooled_assigned,
                              "unassigned_family_n": len(d) - nk})

    for name, rws in [("joint_summary", summaries), ("family_joint_cells", cells),
                      ("joint_support", support), ("genome_system_states", allstates)]:
        pd.DataFrame(rws).to_csv(outdir / (name + ".tsv"), sep="\t", index=False)
    (outdir / "validation.json").write_text(json.dumps(
        {"all_prior_maintenance_calls_match": True, "pairs": 4, "cohorts": 4,
         "CO_activity_claimed": False, "environment_model_fitted": False,
         "thresholds": [2, 3, 5]}, indent=2))
    print(pd.DataFrame(support).query("min_per_joint_state==3").to_string(index=False))


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser(description="Joint-configuration support diagnostic (descriptive only)")
    ap.add_argument("--data-root", type=Path, default=repo / "data/upstream",
                    help="root of staged upstream inputs (default: <repo>/data/upstream)")
    ap.add_argument("--outdir", type=Path, required=True, help="output directory")
    a = ap.parse_args()
    run(a.data_root.resolve(), a.outdir.resolve())


if __name__ == "__main__":
    main()
