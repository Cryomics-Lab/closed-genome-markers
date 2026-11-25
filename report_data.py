"""Small CSVs behind the figures of the cryomics.org report on closed genomes.

  report_data.py <analysis dir> <marker_info dir> <out dir>

<analysis dir> holds markers.py output for both sets (reps.genomes.tsv,
reps.markers.tsv, post232.genomes.tsv, post232.markers.tsv). <marker_info dir>
is GTDB-Tk's align/intermediate_results, for marker names. Writes to <out dir>:
  summary.csv   share of genomes past each error threshold, per set and CheckM2 model
  dist.csv      share of genomes whose error exceeds x, x = 0..20 (%), per set and measure
  markers.csv   per marker: % of genomes missing it / with several copies
  phyla.csv     per phylum (>= 10 genomes in a set): mean errors of both tools
  heat.csv      % of genomes missing each of the most-missed markers, per phylum (reps)
  outliers.csv  the genomes furthest from a perfect score, by each measure
and prints the numbers the text quotes.
"""
import csv, os, sys
from collections import defaultdict, Counter

A, MI, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
SETS = ('reps', 'post232')
info = {}
for s in ('bac120', 'ar53'):
    for r in csv.DictReader(open(f'{MI}/gtdbtk.{s}.marker_info.tsv'), delimiter='\t'):
        info[r['Marker ID']] = (r['Name'], r['Description'])

G = {}
for d in SETS:
    G[d] = []
    for r in csv.DictReader(open(f'{A}/{d}.genomes.tsv'), delimiter='\t'):
        lin = r['lineage'].split(';') if r['lineage'] else []
        n = int(r['n_markers'])
        G[d].append(dict(acc=r['accession'], set=r['marker_set'], n=n,
            miss=int(r['missing']), multi=int(r['multi']),
            miss_pct=100 * int(r['missing']) / n, multi_pct=100 * int(r['multi']) / n,
            deficit=100 - float(r['completeness']), cont=float(r['contamination']),
            model='specific' if r['checkm2_model'].startswith('Neural') else 'general',
            domain='Archaea' if r['marker_set'] == 'ar53' else 'Bacteria',
            phylum=lin[1][3:] if len(lin) > 1 and lin[1] != 'p__' else '',
            species=lin[6][3:] if len(lin) > 6 else (lin[-1][3:] if lin else ''),
            density=float(r['coding_density']), size=int(r['genome_size'])))

def w(name, rows, cols):
    with open(f'{OUT}/{name}', 'w', newline='') as f:
        x = csv.writer(f); x.writerow(cols)
        for r in rows: x.writerow([round(r[c], 2) if isinstance(r[c], float) else r[c] for c in cols])

TESTS = [('compl_lt99', lambda g: g['deficit'] > 1), ('compl_lt95', lambda g: g['deficit'] > 5),
         ('compl_lt90', lambda g: g['deficit'] > 10), ('cont_gt1', lambda g: g['cont'] > 1),
         ('cont_gt2', lambda g: g['cont'] > 2), ('cont_gt5', lambda g: g['cont'] > 5),
         ('checkm2_clean', lambda g: g['deficit'] <= 1 and g['cont'] <= 1),
         ('miss_gt0', lambda g: g['miss'] > 0), ('miss_gt5', lambda g: g['miss'] > 5),
         ('multi_gt0', lambda g: g['multi'] > 0), ('multi_gt5', lambda g: g['multi'] > 5),
         ('markers_clean', lambda g: g['miss'] == 0 and g['multi'] == 0)]
summ = []
for d in SETS:
    for m in ('all', 'specific', 'general'):
        S = [g for g in G[d] if m == 'all' or g['model'] == m]
        row = dict(set=d, model=m, n=len(S))
        for k, t in TESTS: row[k] = 100 * sum(map(t, S)) / len(S)
        summ.append(row)
w('summary.csv', summ, ['set', 'model', 'n'] + [k for k, _ in TESTS])

MEAS = [('deficit', 'CheckM2 completeness below 100'), ('cont', 'CheckM2 contamination'),
        ('miss_pct', 'markers missing'), ('multi_pct', 'markers in several copies')]
dist = []
for d in SETS:
    for k, _ in MEAS:
        v = [g[k] for g in G[d]]
        for i in range(0, 81):
            x = i / 4
            dist.append(dict(set=d, measure=k, x=x, share=100 * sum(y > x for y in v) / len(v)))
w('dist.csv', dist, ['set', 'measure', 'x', 'share'])

mk = []
for d in SETS:
    for r in csv.DictReader(open(f'{A}/{d}.markers.tsv'), delimiter='\t'):
        nm, desc = info.get(r['marker'], ('', ''))
        mk.append(dict(set=d, marker_set=r['marker_set'], marker=r['marker'], name=nm, description=desc,
                       genomes=int(r['genomes']), pct_missing=float(r['pct_missing']), pct_multi=float(r['pct_multi'])))
w('markers.csv', mk, ['set', 'marker_set', 'marker', 'name', 'description', 'genomes', 'pct_missing', 'pct_multi'])

ph = []
for d in SETS:
    by = defaultdict(list)
    for g in G[d]:
        if g['phylum']: by[(g['domain'], g['phylum'])].append(g)
    for (dom, p), S in by.items():
        if len(S) < 10: continue
        mean = lambda k: sum(g[k] for g in S) / len(S)
        ph.append(dict(set=d, domain=dom, phylum=p, n=len(S), miss_pct=mean('miss_pct'), multi_pct=mean('multi_pct'),
                       deficit=mean('deficit'), cont=mean('cont'),
                       general=100 * sum(g['model'] == 'general' for g in S) / len(S),
                       miss_gt5=100 * sum(g['miss'] > 5 for g in S) / len(S),
                       cont_gt5=100 * sum(g['cont'] > 5 for g in S) / len(S)))
ph.sort(key=lambda r: (r['set'], -r['n']))
w('phyla.csv', ph, ['set', 'domain', 'phylum', 'n', 'miss_pct', 'multi_pct', 'deficit', 'cont', 'general', 'miss_gt5', 'cont_gt5'])

# heat map: reps only (post232 has too few genomes outside a handful of phyla)
miss_lists = {}
for r in csv.DictReader(open(f'{A}/reps.genomes.tsv'), delimiter='\t'):
    miss_lists[r['accession']] = set(filter(None, r['missing_list'].split(',')))
heat = []
for ms, dom, nmark, nmin in (('bac120', 'Bacteria', 30, 20), ('ar53', 'Archaea', 20, 10)):
    top = [r for r in mk if r['set'] == 'reps' and r['marker_set'] == ms]
    top = [r['marker'] for r in sorted(top, key=lambda r: -r['pct_missing'])[:nmark]]
    by = defaultdict(list)
    for g in G['reps']:
        if g['domain'] == dom and g['phylum']: by[g['phylum']].append(g)
    for p, S in by.items():
        if len(S) < nmin: continue
        for m in top:
            heat.append(dict(marker_set=ms, phylum=p, n=len(S), marker=m, name=info[m][0],
                             pct_missing=100 * sum(m in miss_lists[g['acc']] for g in S) / len(S)))
w('heat.csv', heat, ['marker_set', 'phylum', 'n', 'marker', 'name', 'pct_missing'])

out = []
for d in SETS:
    for k, _ in (('cont', 0), ('deficit', 0), ('miss_pct', 0)):
        for g in sorted(G[d], key=lambda g: -g[k])[:10]:
            out.append(dict(set=d, by=k, accession=g['acc'], phylum=g['phylum'], species=g['species'],
                            completeness=100 - g['deficit'], contamination=g['cont'], model=g['model'],
                            missing=g['miss'], multi=g['multi'], n_markers=g['n'], coding_density=g['density']))
w('outliers.csv', out, ['set', 'by', 'accession', 'phylum', 'species', 'completeness', 'contamination', 'model',
                        'missing', 'multi', 'n_markers', 'coding_density'])

# numbers quoted in the text
def spearman(x, y):
    def rank(v):
        o = sorted(range(len(v)), key=v.__getitem__); r = [0.0] * len(v); i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]: j += 1
            for k in range(i, j + 1): r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    a, b = rank(x), rank(y); n = len(a); ma, mb = sum(a) / n, sum(b) / n
    cov = sum((p - ma) * (q - mb) for p, q in zip(a, b))
    return cov / (sum((p - ma) ** 2 for p in a) * sum((q - mb) ** 2 for q in b)) ** .5
for d in SETS:
    S = G[d]
    print(f"== {d}: {len(S)} genomes, {sum(g['domain'] == 'Archaea' for g in S)} archaeal, "
          f"{sum(not g['phylum'] for g in S)} without a phylum, {Counter(g['model'] for g in S)}")
    print(f"   spearman(markers missing %, completeness deficit) = {spearman([g['miss_pct'] for g in S], [g['deficit'] for g in S]):.2f}")
    print(f"   spearman(markers multi %, contamination) = {spearman([g['multi_pct'] for g in S], [g['cont'] for g in S]):.2f}")
for r in summ: print('  ', {k: (round(v, 1) if isinstance(v, float) else v) for k, v in r.items()})
