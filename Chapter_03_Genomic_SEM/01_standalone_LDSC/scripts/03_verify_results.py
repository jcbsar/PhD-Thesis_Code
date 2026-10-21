#!/usr/bin/env python3
"""Validate the standalone LDSC results reported in Figure 3.1.

Checks SNP-based heritability, genetic correlations, confidence intervals,
and BH-FDR-adjusted p-values using summary-level results.

Optionally compares estimates against the original LDSC execution logs.
Uses only Python's standard library and does not rerun LDSC.
"""
from __future__ import annotations
import argparse
import csv
import math
import re
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / 'data'
OUT = BASE / 'results'
EPI = {'CAE', 'GTCS', 'JAE', 'JME', 'F-HS', 'F-Neg', 'F-Oth'}
BRAIN = {'ACC', 'AMY', 'BRS', 'CAU', 'PAL', 'PUT', 'THA'}
ORDER = ['ACC', 'AMY', 'BRS', 'CAU', 'PAL', 'PUT', 'THA',
         'CAE', 'GTCS', 'JAE', 'JME', 'F-HS', 'F-Neg', 'F-Oth']
LOG_NAMES = {'CAE': 'CAE', 'GTCS': 'GTCS', 'JAE': 'JAE', 'JME': 'JME',
             'FE_HS': 'F-HS', 'FE_Neg': 'F-Neg', 'FE_Oth': 'F-Oth',
             'accumbens': 'ACC', 'amygdala': 'AMY', 'brainstem': 'BRS',
             'caudate': 'CAU', 'pallidum': 'PAL', 'putamen': 'PUT',
             'thalamus': 'THA'}


def rows(path: Path, delimiter=','):
    with path.open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def pair_key(a, b):
    return tuple(sorted((a, b)))


def check(condition, reason):
    if not condition:
        raise AssertionError(reason)


def approx(a, b, tol=0.00011):
    return math.isclose(float(a), float(b), abs_tol=tol, rel_tol=0)


def fdr_bh(pvalues):
    """Benjamini–Hochberg q values in original order (no dependencies)."""
    n = len(pvalues)
    indices = sorted(range(n), key=lambda i: pvalues[i])
    out = [float('nan')] * n
    minimum = 1.0
    for rank_position in range(n - 1, -1, -1):
        idx = indices[rank_position]
        minimum = min(minimum, pvalues[idx] * n / (rank_position + 1), 1.0)
        out[idx] = minimum
    return out


def check_heritability(log_dir=None):
    df = rows(DATA / 'figure_3_1A_heritability.csv')
    check(len(df) == 16, 'Expected 16 heritability entries: 9 epilepsy and 7 volumes')
    name_map = {r['Phenotype']: r for r in df}
    check(len(name_map) == len(df), 'Duplicate heritability phenotype')
    check(set(name_map) == EPI | BRAIN | {'GGE', 'FE'}, 'Unexpected heritability phenotypes')
    for r in df:
        x, se = float(r['h2']), float(r['SE'])
        check(approx(float(r['CI_lower']), x - 1.96 * se, 0.0000001), 'Incorrect lower CI for ' + r['Phenotype'])
        check(approx(float(r['CI_upper']), x + 1.96 * se, 0.0000001), 'Incorrect upper CI for ' + r['Phenotype'])
    check(approx(name_map['GGE']['h2'], 0.2317, 1e-8), 'Unexpected GGE value')
    check(approx(name_map['FE']['h2'], 0.0287, 1e-8), 'Unexpected FE value')
    if log_dir:
        count = 0
        for stem, display in LOG_NAMES.items():
            # Standalone --h2 logs used FE_HS.log etc and brain structure names.
            logfile = log_dir / (stem + '.log')
            check(logfile.is_file(), 'Missing optional h2 log: ' + str(logfile))
            content = logfile.read_text(encoding='utf-8')
            match = re.search(r'Total (?:Liability|Observed) scale h2:\s*([0-9.eE+-]+)\s*\(([0-9.eE+-]+)\)', content)
            check(match is not None, 'No h2/SE pair in ' + str(logfile))
            check(approx(match[1], name_map[display]['h2']), 'h2 disagreement for ' + stem)
            check(approx(match[2], name_map[display]['SE']), 'SE disagreement for ' + stem)
            count += 1
        print(f'Validated {count} archived univariate LDSC logs against Figure 3.1A.')
    else:
        print('Validated Figure 3.1A exported values and confidence intervals; original logs not supplied to this run.')
    print('Broad GGE and FE are archive-derived values; original --h2 logs were not supplied.')
    return df


def check_correlations(log_dir=None):
    rg = rows(DATA / 'figure_3_1B_rg_summary.tsv', '\t')
    check(len(rg) == 91, 'Expected exactly 91 unique correlations among 14 phenotypes')
    by_pair = {}
    for r in rg:
        a, b = r['trait1'], r['trait2']
        check(a in ORDER and b in ORDER and a != b, 'Invalid pair: '+a+', '+b)
        key = pair_key(a, b)
        check(key not in by_pair, 'Duplicate pair '+str(key))
        by_pair[key] = r
    check(len(by_pair) == 14 * 13 // 2, 'Incomplete genetic-correlation table')

    matrix = list(csv.reader((DATA / 'figure_3_1B_rg_matrix.csv').open(newline='', encoding='utf-8-sig')))
    check(len(matrix) == 15 and len(matrix[0]) == 15, 'Expected 14x14 rg matrix')
    check(matrix[0][1:] == ORDER, 'Unexpected column order in rg matrix')
    for i, row in enumerate(matrix[1:]):
        check(row[0] == ORDER[i] and len(row) == 15, 'Unexpected row in rg matrix')
        for j, val in enumerate(row[1:]):
            if i == j:
                check(approx(val, 1, 0), 'Matrix diagonal is not 1')
            else:
                check(approx(val, by_pair[pair_key(ORDER[i], ORDER[j])]['rg'], 0.0000001),
                      'Matrix versus long-format rg discrepancy')

    cross = [r for r in rg if (r['trait1'] in EPI and r['trait2'] in BRAIN) or
                                  (r['trait2'] in EPI and r['trait1'] in BRAIN)]
    check(len(cross) == 49, 'Expected 49 cross-domain comparisons')
    qvals = fdr_bh([float(r['p']) for r in cross])
    nominal = [(r['trait1'], r['trait2']) for r in cross if float(r['p']) < .05]
    expected_nominal = {('GTCS', 'AMY'), ('JAE', 'CAU'), ('JME', 'BRS'), ('JME', 'PAL')}
    check(set(nominal) == expected_nominal, 'Nominal associations differ from thesis text')
    check(all(q >= .05 for q in qvals), 'Unexpected FDR-significant association')
    check(all(approx(q, 0.521850375, 0.00001) for r, q in zip(cross, qvals)
              if (r['trait1'], r['trait2']) in expected_nominal), 'Unexpected q for nominal pairs')

    OUT.mkdir(exist_ok=True)
    destination = OUT / 'figure_3_1_cross_domain_fdr.csv'
    with destination.open('w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['epilepsy_phenotype', 'brain_volume', 'rg', 'se', 'z', 'p_unadjusted', 'q_BH_FDR'])
        for r, q in zip(cross, qvals):
            ep = r['trait1'] if r['trait1'] in EPI else r['trait2']
            bv = r['trait2'] if r['trait2'] in BRAIN else r['trait1']
            writer.writerow([ep, bv, r['rg'], r['se'], r['z'], r['p'], format(q, '.12g')])

    if log_dir:
        # Original standalone rg logs contain a "Summary of Genetic Correlation Results"
        # table. Read only its pairwise lines, not the full path-identifying header.
        observed = {}
        for logfile in sorted(log_dir.glob('*_rg.log')):
            lines = logfile.read_text(encoding='utf-8').splitlines()
            try:
                marker = next(i for i, l in enumerate(lines) if l.startswith('Summary of Genetic Correlation Results'))
            except StopIteration:
                raise AssertionError('Missing summary table in ' + str(logfile))
            for line in lines[marker+2:]:
                parts = line.split()
                if len(parts) < 6 or not parts[0].endswith('.sumstats.gz') or not parts[1].endswith('.sumstats.gz'):
                    continue
                stems = [Path(path).name.replace('.hm3.sumstats.gz', '') for path in parts[:2]]
                a, b = [next((display for original, display in LOG_NAMES.items()
                              if stem == original or stem.startswith(original+'_cau_') or stem.startswith(original+'_21_') or
                              (original == 'FE_HS' and stem.startswith('focal_HS_')) or
                              (original == 'FE_Neg' and stem.startswith('focal_lesion_negative_')) or
                              (original == 'FE_Oth' and stem.startswith('focal_lesion_other_')) or
                              (original == 'GTCS' and stem.startswith('GTCS_cau_'))), None)
                        for stem in stems]
                # Brain GWAS filenames include a date stamp; map by name prefix.
                check(a is not None and b is not None, 'Unknown original name(s): '+str(stems))
                key = pair_key(a, b)
                check(key not in observed, 'Duplicate archived rg: '+str(key))
                observed[key] = [float(v) for v in parts[2:6]]
        check(len(observed) == 91, f'Expected 91 log estimates, found {len(observed)}')
        for key, (value, se, z, pval) in observed.items():
            ref = by_pair[key]
            for label, actual in zip(['rg', 'se', 'z', 'p'], [value, se, z, pval]):
                check(approx(actual, ref[label], 0.0000001), f'{key}: mismatched {label}')
        print('Validated 91 pairwise estimates (rg/SE/Z/p) against 13 archived LDSC logs.')
    else:
        print('Validated 91 correlation values against saved matrix; original logs not supplied to this run.')
    print('BH FDR: 49 cross-domain tests, four nominal associations, zero q < 0.05.')
    print('Wrote derived results to ' + str(destination))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--h2-logs', type=Path, help='private folder with original 14 standalone --h2 .log files')
    parser.add_argument('--rg-logs', type=Path, help='private folder with original 13 standalone --rg .log files')
    args = parser.parse_args()
    check_heritability(args.h2_logs)
    check_correlations(args.rg_logs)
    print('PASS: archived summary-level outputs are internally consistent.')


if __name__ == '__main__':
    main()
