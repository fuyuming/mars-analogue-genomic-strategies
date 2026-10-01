"""Same-raw-output parser sensitivity, not an ecological CaCo replication.

Input: completed HMMER3 hmmsearch domtblout and frozen 16-genome manifest.
Output: per-genome annotation counts, raw accepted domains and set differences.
No pairwise competition, environmental association or imputed absence is computed.
"""
from pathlib import Path
import ast
import collections
import hashlib
from itertools import combinations
import json
import os
import pandas as pd

H = Path(__file__).resolve().parent
O = H / 'host/pilot_comparison'


def parser(version, coverage_only=False):
    src = H / 'upstream' / version / 'CaCo.py'
    name = 'process_hmmsearch_output' if version == 'archived' else 'parse_hmmsearch_output'
    f = next(n for n in ast.parse(src.read_text()).body
             if isinstance(n, ast.FunctionDef) and n.name == name)
    if coverage_only:
        # Change only qlen=int(fields[2]) to qlen=int(fields[5]).
        assignments = [n for n in ast.walk(f) if isinstance(n, ast.Assign)
                       and any(isinstance(t, ast.Name) and t.id == 'qlen' for t in n.targets)]
        assert len(assignments) == 1
        assignments[0].value.args[0].slice = ast.Constant(value=5)
    module = ast.fix_missing_locations(ast.Module(body=[f], type_ignores=[]))
    ns = {'pd': pd, 'os': os}
    exec(compile(module, str(src), 'exec'), ns)
    return ns[name]


def base_family(model):
    return model.split('.')[0].split('_')[0]


def substrates(families, mapping, exclude_cbm=False):
    return {s.strip() for family in families if family in mapping
            and not (exclude_cbm and family.startswith('CBM'))
            for s in mapping[family].split(',') if s.strip()}


def main():
    O.mkdir(exist_ok=True)
    raw = H / 'host/pilot16.domtbl'
    manifest = pd.read_csv(H / 'pilot16_manifest.tsv', sep='\t')
    names = set(manifest.Name)
    assert len(names) == 16
    assert '# [ok]' in raw.read_text()[-1000:]
    receipt = json.loads((H / 'host/hmmsearch_complete.json').read_text())
    assert hashlib.sha256(raw.read_bytes()).hexdigest() == receipt['raw_sha256']
    mapping = json.loads((H / 'upstream/archived/data__substrate_key.json').read_text())
    lines = {n: [] for n in names}
    records = []
    for line in raw.open():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split()
        assert len(f) >= 22
        genome = f[0].split('|')[0]
        assert genome in names, genome
        lines[genome].append(line)
        qlen, hs, he = int(f[5]), int(f[15]), int(f[16])
        assert 1 <= hs <= he <= qlen
        records.append(dict(genome=genome, protein=f[0], model=f[3], family=base_family(f[3]),
                            protein_length=int(f[2]), model_length=qlen,
                            seq_e=float(f[6]), seq_score=float(f[7]), domain_e=float(f[12]),
                            domain_score=float(f[13]), hmm_from=hs, hmm_to=he,
                            ali_from=int(f[17]), ali_to=int(f[18]),
                            model_coverage=(he-hs+1)/qlen,
                            upstream_coverage=min(1., (he-hs+1)/int(f[2]))))
    domains = pd.DataFrame(records)
    assert len(domains) > 0
    domains['seq_gate'] = ((domains.seq_e <= 1e-15) & (domains.seq_score >= 25)
                           & (domains.model_coverage >= .35))
    domains['domain_gate'] = ((domains.domain_e <= 1e-15) & (domains.domain_score >= 25)
                              & (domains.model_coverage >= .35))
    domains.to_csv(O / 'raw_domains_with_gates.tsv', sep='\t', index=False)
    conflicts = []
    for (genome, protein), group in domains[domains.domain_gate].groupby(['genome', 'protein']):
        for a, b in combinations(group.itertuples(), 2):
            if a.family == b.family:
                continue
            overlap = max(0, min(a.ali_to, b.ali_to) - max(a.ali_from, b.ali_from) + 1)
            fraction = overlap / min(a.ali_to-a.ali_from+1, b.ali_to-b.ali_from+1)
            if fraction >= .5:
                conflicts.append(dict(genome=genome, protein=protein, family_a=a.family,
                                      family_b=b.family, overlap_fraction_of_shorter=fraction))
    pd.DataFrame(conflicts, columns=['genome','protein','family_a','family_b','overlap_fraction_of_shorter']).to_csv(
        O / 'overlapping_family_candidates.tsv', sep='\t', index=False)
    policies = [('archived', parser('archived')), ('current', parser('current')),
                ('coverage_only', parser('current', True))]
    accepted = collections.defaultdict(dict)
    summaries = []
    diffs = []
    for genome in sorted(names):
        sub = domains[domains.genome == genome]
        for label, func in policies:
            directory = O / label
            directory.mkdir(exist_ok=True)
            p = directory / (genome + '.domtbl')
            p.write_text('# Per-genome subset of a single joint search, preserving HMMER order\n'
                         + ''.join(lines[genome]))
            parsed = Path(str(p) + '.parsed')
            if parsed.exists():
                parsed.unlink()  # reproducible derived file, never a raw search result
            func(str(p), str(directory))
            x = pd.read_csv(parsed, sep='\t')
            pairs = {(r.Query_Domain, base_family(r.Hit_Name)) for r in x.itertuples()}
            accepted[label][genome] = pairs
        for label, gate in [('all_seq_gate', 'seq_gate'), ('all_domain_gate', 'domain_gate')]:
            accepted[label][genome] = set(zip(sub.loc[sub[gate], 'protein'], sub.loc[sub[gate], 'family']))
        for label in accepted:
            pairs = accepted[label][genome]
            fams = {f for _, f in pairs}
            prots = {p for p, _ in pairs}
            per_protein = collections.Counter(p for p, _ in pairs)
            summaries.append(dict(genome=genome, cohort=manifest.set_index('Name').loc[genome, 'cohort'],
                                  policy=label, proteins=len(prots), protein_family_pairs=len(pairs),
                                  families=len(fams), mapped_families=len(fams & mapping.keys()),
                                  substrate_categories=len(substrates(fams, mapping)),
                                  substrates_without_cbm_only=len(substrates(fams, mapping, True)),
                                  proteins_multiple_families=sum(n > 1 for n in per_protein.values())))
        cur, cov, allseq = (accepted[k][genome] for k in ['current', 'coverage_only', 'all_seq_gate'])
        diffs.append(dict(genome=genome, current_pairs=len(cur), coverage_only_pairs=len(cov),
                          gained_by_coverage=len(cov-cur), lost_by_coverage=len(cur-cov),
                          allseq_pairs=len(allseq), gained_by_preserving_families=len(allseq-cov),
                          lost_after_preserving_families=len(cov-allseq),
                          current_substrates=len(substrates({f for _, f in cur}, mapping)),
                          allseq_substrates=len(substrates({f for _, f in allseq}, mapping))))
    summary = pd.DataFrame(summaries)
    summary.to_csv(O / 'per_genome_policy_counts.tsv', sep='\t', index=False)
    pd.DataFrame(diffs).to_csv(O / 'per_genome_policy_differences.tsv', sep='\t', index=False)
    # Store exact pairs, so changes in breadth can be traced to raw domain rows.
    pd.DataFrame([dict(policy=k, genome=g, protein=p, family=f)
                  for k, v in accepted.items() for g, pairs in v.items() for p, f in sorted(pairs)]).to_csv(
        O / 'accepted_protein_family_pairs.tsv', sep='\t', index=False)
    evidence = []
    for genome, pairs in accepted['all_domain_gate'].items():
        counts = collections.Counter(f for _, f in pairs)
        for family, count in sorted(counts.items()):
            for substrate in sorted(substrates({family}, mapping)):
                evidence.append(dict(genome=genome, family=family, substrate=substrate,
                                     proteins=count, binding_module=family.startswith('CBM')))
    evidence = pd.DataFrame(evidence).sort_values(['genome','family','substrate'])
    evidence.to_csv(O / 'substrate_family_evidence.tsv', sep='\t', index=False)
    support = evidence[~evidence.binding_module].groupby(['genome','substrate']).agg(
        supporting_families=('family','nunique'), protein_assignments=('proteins','sum')).reset_index()
    support.to_csv(O / 'substrate_support_counts.tsv', sep='\t', index=False)
    all_pairs = [(p, f) for pairs in accepted['all_domain_gate'].values() for p, f in pairs]
    observed_families = {f for _, f in all_pairs}
    mapping_report = dict(domain_gate_protein_family_pairs=len(all_pairs),
                          mapped_pairs=sum(f in mapping for _, f in all_pairs),
                          genome_substrate_assignments_excluding_cbm_only=len(support),
                          supported_by_one_family=int(support.supporting_families.eq(1).sum()),
                          one_family_percent=round(100*support.supporting_families.eq(1).mean(), 2),
                          unique_substrate_labels=sorted(evidence.substrate.unique().tolist()),
                          observed_families_mapping_to_multiple_substrates=sum(',' in mapping.get(f,'') for f in observed_families),
                          observed_mapped_families=sum(f in mapping for f in observed_families),
                          boundary='Family support counts quantify mapping redundancy, not independent functional validation or confirmed uptake.')
    (O / 'substrate_evidence_summary.json').write_text(json.dumps(mapping_report, indent=2))
    coverage_false_negative = (domains.seq_e <= 1e-15) & (domains.seq_score >= 25) & (domains.model_coverage >= .35) & (domains.upstream_coverage < .35)
    coverage_false_positive = (domains.seq_e <= 1e-15) & (domains.seq_score >= 25) & (domains.model_coverage < .35) & (domains.upstream_coverage >= .35)
    aggregate = summary.groupby('policy')[['proteins','protein_family_pairs','families','substrate_categories','substrates_without_cbm_only','proteins_multiple_families']].sum().astype(int).to_dict('index')
    outcome = dict(genomes=16, raw_domain_rows=len(domains), per_policy_sums=aggregate,
                   proteins_with_overlapping_family_candidates=len({r['protein'] for r in conflicts}),
                   coverage_discordant_domain_rows=dict(upstream_reject_model_accept=int(coverage_false_negative.sum()),
                                                        upstream_accept_model_reject=int(coverage_false_positive.sum())),
                   genomes_changed_current_vs_allseq=int(sum(accepted['current'][n] != accepted['all_seq_gate'][n] for n in names)),
                   genomes_with_changed_substrate_sets=int(sum(substrates({f for _, f in accepted['current'][n]}, mapping) != substrates({f for _, f in accepted['all_seq_gate'][n]}, mapping) for n in names)),
                   raw_sha256=receipt['raw_sha256'],
                   limitations=['Annotation method pilot only; no ecological hypothesis test.',
                                'Substrate mappings are candidate capabilities, not verified consumption.',
                                'All-family policies retain possible overlapping/ambiguous assignments.',
                                'Aggregated families/substrates sum per-genome counts, not unique totals.',
                                'Joint target database search is not a full per-genome CaCo reproduction.'])
    (O / 'comparison_summary.json').write_text(json.dumps(outcome, indent=2))
    print(json.dumps(outcome, indent=2))


if __name__ == '__main__':
    main()
