"""Does GTDB-Tk's alignment drop markers found in differing copies?

  msa_check.py <run out dir>

For every bacterial genome in the alignment, compares the filled fraction of
its row in gtdbtk.bac120.user_msa.fasta.gz with the fraction of bac120 markers
found once or in identical copies. A Pearson r near 1 means a marker in
differing copies is as absent from the alignment as a missing one.
"""
import csv, gzip, sys
out = sys.argv[1]
summ = {}
for r in csv.DictReader(open(f'{out}/gtdbtk/identify/gtdbtk.bac120.markers_summary.tsv'), delimiter='\t'):
    summ[r.get('name') or r['Name']] = r
seqs, cur = {}, None
for l in gzip.open(f'{out}/gtdbtk/align/gtdbtk.bac120.user_msa.fasta.gz', 'rt'):
    if l[0] == '>': cur = l[1:].strip(); seqs[cur] = []
    else: seqs[cur].append(l.strip())
xs, ys = [], []
for k, v in seqs.items():
    s = ''.join(v); r = summ[k]
    n = sum(int(r[c]) for c in ('number_unique_genes', 'number_multiple_genes',
                                'number_multiple_unique_genes', 'number_missing_genes'))
    xs.append((int(r['number_unique_genes']) + int(r['number_multiple_unique_genes'])) / n)
    ys.append(sum(c not in '-.' for c in s) / len(s))
mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
r = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / (
    sum((a - mx) ** 2 for a in xs) * sum((b - my) ** 2 for b in ys)) ** .5
print(f'{len(xs)} genomes, alignment length {len(s)}: Pearson r(usable markers, filled alignment) = {r:.3f}')
