"""Give genomes outside GTDB a GTDB lineage from their NCBI name.

  ncbi_to_gtdb_lineage.py <bac120_metadata.tsv.gz> <ar53_metadata.tsv.gz> <picks.tsv> > lineage.tsv

<picks.tsv> has a header and columns accession, species_taxid, organism (as
written by select_post232.py). Every GTDB genome pairs an NCBI lineage with a
GTDB lineage. A genome's NCBI species, or failing that the first word of its
name as an NCBI genus, family, order, class or phylum ("Woeseiaceae bacterium"),
takes the GTDB lineage its GTDB genomes share, down to the deepest rank on which
at least 90% of them agree. Output: accession, lineage, matched_on
(species|g|f|o|c|p|none).
Unmatched genomes get no lineage and drop out of taxonomy-wise summaries.
"""
import csv, gzip, sys
from collections import Counter, defaultdict

AGREE = 0.9
RANKS = 'gfocp'
by_sp, by_rank = defaultdict(Counter), {k: defaultdict(Counter) for k in RANKS}
for f in sys.argv[1:3]:
    for r in csv.DictReader(gzip.open(f, 'rt'), delimiter='\t'):
        nt = dict(x.split('__', 1) for x in r['ncbi_taxonomy'].split(';') if '__' in x)
        g = tuple(r['gtdb_taxonomy'].split(';'))
        if nt.get('s'): by_sp[nt['s']][g] += 1
        for k in RANKS:
            if nt.get(k): by_rank[k][nt[k]][g] += 1

def consensus(c):
    """Deepest GTDB lineage prefix shared by at least AGREE of the genomes in c."""
    n, best = sum(c.values()), ()
    for depth in range(1, 8):
        pre = Counter()
        for g, k in c.items(): pre[g[:depth]] += k
        top, k = pre.most_common(1)[0]
        if k / n < AGREE: break
        best = top
    return ';'.join(best)

def binomial(name):
    w = [x for x in name.split() if x not in ('uncultured', 'Candidatus')]
    return ' '.join(w[:2]), (w[0] if w else '')

print('accession\tlineage\tmatched_on')
for r in csv.DictReader(open(sys.argv[3]), delimiter='\t'):
    sp, first = binomial(r['organism'])
    for key, table, how in [(sp, by_sp, 'species')] + [(first, by_rank[k], k) for k in RANKS]:
        if key in table and (lin := consensus(table[key])):
            print(f"{r['accession']}\t{lin}\t{how}"); break
    else:
        print(f"{r['accession']}\t\tnone")
