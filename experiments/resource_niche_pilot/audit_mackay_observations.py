"""Audit missing measurements without converting dashes to biological absence."""
from pathlib import Path
import argparse
import json
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--coverage', type=Path, required=True)
    p.add_argument('--parents', type=Path, required=True)
    p.add_argument('--quality', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    cov = pd.read_csv(a.coverage, sep='\t', keep_default_na=False)
    parents = pd.read_csv(a.parents, sep='\t')
    quality = pd.read_csv(a.quality, sep='\t')
    assert not cov.duplicated(['MAG','metagenome_id']).any()
    assert len(cov) == 451 * 16 and cov.MAG.nunique() == 451
    assert not parents.isolate.duplicated().any()
    assert not quality.Name.duplicated().any()
    raw = pd.to_numeric(cov.metagenome_percent, errors='coerce')
    assert set(cov.loc[raw.isna(), 'metagenome_percent']) == {'-'}
    cov['numeric'] = raw.notna()
    cov['explicit_zero'] = raw.eq(0)
    counts = cov.groupby('MAG').agg(n_numeric=('numeric','sum'), n_explicit_zero=('explicit_zero','sum'))
    counts['n_dash'] = 16 - counts.n_numeric
    parents['Name'] = 'NEE56_mackay__' + parents.assembly_accession
    quality['primary_recalc'] = (quality.Completeness >= 50) & (quality.Contamination < 10)
    quality['strict_recalc'] = (quality.Completeness >= 90) & (quality.Contamination < 5)
    joined = parents.merge(quality[['Name','Completeness','Contamination','primary_recalc','strict_recalc']],
                           on='Name', how='left', validate='one_to_one', indicator=True)
    assert joined['_merge'].eq('both').all()
    joined = joined.drop(columns='_merge').merge(counts, left_on='isolate', right_index=True, validate='one_to_one')
    assert len(joined) == 451
    assert (joined.n_numeric == joined.n_assembly_parents).all()
    joined.to_csv(a.out / 'observation_audit_reproduced.tsv', sep='\t', index=False)
    summary = dict(MAGs=451, samples=16, numeric_entries=int(cov.numeric.sum()),
                   missing_dash_entries=int((~cov.numeric).sum()), explicit_zero_entries=int(cov.explicit_zero.sum()),
                   all_numeric_per_MAG_equals_assembly_parent_count=True,
                   by_assembly=joined.groupby('n_assembly_parents').agg(
                       total=('Name','size'), primary=('primary_recalc','sum'), strict=('strict_recalc','sum')
                   ).reset_index().to_dict('records'))
    (a.out / 'observation_summary_reproduced.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
