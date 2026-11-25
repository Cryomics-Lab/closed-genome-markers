"""Which genomes are closed: NCBI's topology for every sequence.

  topology.py <workdir>

Reads <workdir>/seqs.tsv (accession, sequence id, header), asks NCBI
E-utilities esummary (db=nuccore) for each sequence's topology in batches of
200, caching the replies in <workdir>/esum/, and writes
<workdir>/genome_topology.tsv: accession, n_seqs, n_circular, largest_seq,
largest_len, largest_topology, all_circular. A genome counts as closed when its
largest sequence is circular. Sequences esummary does not return get topology NA.
"""
import glob, json, os, sys, time, urllib.parse, urllib.request
from collections import Counter, defaultdict

W = sys.argv[1]
os.makedirs(f'{W}/esum', exist_ok=True)
g = defaultdict(list)
for l in open(f'{W}/seqs.tsv'):
    f = l.rstrip('\n').split('\t'); g[f[0]].append(f[1])
ids = sorted({s for ss in g.values() for s in ss})
for i in range(0, len(ids), 200):
    out = f'{W}/esum/chunk_{i // 200:04d}.json'
    if os.path.exists(out) and os.path.getsize(out): continue
    body = urllib.parse.urlencode({'db': 'nuccore', 'retmode': 'json', 'tool': 'gtdb_topology',
                                   'id': ','.join(ids[i:i + 200])}).encode()
    for t in range(5):
        try:
            d = urllib.request.urlopen('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi',
                                       data=body, timeout=120).read(); break
        except Exception as e:
            print('retry', i, e, file=sys.stderr); time.sleep(5)
    else: sys.exit(f'esummary batch {i} failed')
    open(out, 'wb').write(d); time.sleep(0.4)          # NCBI allows 3 requests/s without a key
top = {}
for f in glob.glob(f'{W}/esum/*.json'):
    d = json.load(open(f))['result']
    for u in d['uids']:
        r = d[u]; top[r['accessionversion']] = (r.get('topology', ''), int(r.get('slen', 0)))
cls = Counter()
with open(f'{W}/genome_topology.tsv', 'w') as o:
    o.write('accession\tn_seqs\tn_circular\tlargest_seq\tlargest_len\tlargest_topology\tall_circular\n')
    for a, ss in g.items():
        rec = [(s,) + top.get(s, ('NA', 0)) for s in ss]
        big = max(rec, key=lambda r: r[2]); nc = sum(r[1] == 'circular' for r in rec)
        o.write(f'{a}\t{len(rec)}\t{nc}\t{big[0]}\t{big[2]}\t{big[1]}\t{nc == len(rec)}\n')
        cls['largest sequence ' + (big[1] or 'blank')] += 1
print(f'{len(g)} genomes, {len(top)} of {len(ids)} sequences with a topology record')
for k, v in sorted(cls.items()): print(f'{v:7d}  {k}')
