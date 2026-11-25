"""Marker-gene and CheckM2 errors on closed genomes.

  markers.py <run out dir> <taxonomy tsv: accession, lineage> <prefix>

Each genome is scored against its own domain's marker set (bac120 or ar53),
taken from the lineage's d__ field; a genome with no lineage takes the set in
which fewer markers are missing. A closed genome should have every marker
exactly once, 100% completeness and 0% contamination; anything else is an
error. Writes <prefix>.genomes.tsv (per genome) and <prefix>.markers.tsv
(per marker: how often missing or multi-copy), and prints summaries.
"""
import csv, sys
from collections import defaultdict, Counter
out, taxf, pre = sys.argv[1:4]
tax = {}
for l in open(taxf):
    f = l.rstrip('\n').split('\t')
    if f[0] != 'accession': tax[f[0]] = f[1].split(';')
def acc(n): return n.removeprefix('q_').removesuffix('.fna').removesuffix('.fna.gz')
ck = {acc(r['Name']): r for r in csv.DictReader(open(f'{out}/checkm2/quality_report.tsv'), delimiter='\t')}
mk = {}
for s in ('bac120', 'ar53'):
    for r in csv.DictReader(open(f'{out}/gtdbtk/identify/gtdbtk.{s}.markers_summary.tsv'), delimiter='\t'):
        mk[(acc(r.get('name') or r['Name']), s)] = r
genomes = []
per_marker = defaultdict(Counter)
for a, c in ck.items():
    lin = [x for x in tax.get(a, []) if x]
    if lin: s = 'ar53' if lin[0] == 'd__Archaea' else 'bac120'
    else:
        if (a, 'bac120') not in mk: continue
        s = min(('bac120', 'ar53'), key=lambda s: int(mk[(a, s)]['number_missing_genes']) / (
            sum(int(mk[(a, s)][k]) for k in ('number_unique_genes', 'number_multiple_genes',
                'number_multiple_unique_genes', 'number_missing_genes'))))
    r = mk.get((a, s))
    if r is None: continue
    miss = [m for m in r['list_missing_genes'].split(',') if m]
    mult = [m for m in (r['list_multiple_genes'] + ',' + r['list_multiple_unique_genes']).split(',') if m]
    uniq = [m for m in r['list_unique_genes'].split(',') if m]
    for m in uniq: per_marker[(s, m)]['ok'] += 1
    for m in miss: per_marker[(s, m)]['missing'] += 1
    for m in mult: per_marker[(s, m)]['multi'] += 1
    genomes.append(dict(accession=a, marker_set=s, n_markers=len(uniq) + len(miss) + len(mult),
        missing=len(miss), multi=len(mult), completeness=float(c['Completeness']),
        contamination=float(c['Contamination']), checkm2_model=c['Completeness_Model_Used'].split(' (')[0],
        genome_size=c['Genome_Size'], coding_density=c['Coding_Density'],
        missing_list=','.join(miss), multi_list=','.join(mult), lineage=';'.join(lin)))
with open(f'{pre}.genomes.tsv', 'w') as f:
    w = csv.DictWriter(f, genomes[0].keys(), delimiter='\t'); w.writeheader(); w.writerows(genomes)
with open(f'{pre}.markers.tsv', 'w') as f:
    f.write('marker_set\tmarker\tgenomes\tmissing\tmulti\tpct_missing\tpct_multi\n')
    for (s, m), c in sorted(per_marker.items(), key=lambda x: -(x[1]['missing'] + x[1]['multi'])):
        n = sum(c.values()); f.write(f"{s}\t{m}\t{n}\t{c['missing']}\t{c['multi']}\t{100*c['missing']/n:.2f}\t{100*c['multi']/n:.2f}\n")
n = len(genomes)
print(f'{n} genomes scored ({sum(g["marker_set"]=="ar53" for g in genomes)} archaeal)')
for k in (0, 1, 3, 5, 10):
    print(f'  missing >{k} markers: {sum(g["missing"]>k for g in genomes):5d}   multi-copy >{k}: {sum(g["multi"]>k for g in genomes):5d}')
print(f'  all markers single-copy: {sum(g["missing"]==0 and g["multi"]==0 for g in genomes)}')
