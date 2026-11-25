"""Complete genomes released to GenBank after a GTDB release's NCBI cutoff,
one per species, downloaded as FASTA.

  select_post232.py <workdir> [first release date, default 09/05/2025]

R232 took NCBI genomes up to 2025-09-04. Writes to <workdir>:
  reports.jsonl                 NCBI Datasets reports, Bacteria (2) and Archaea (2157)
  post232_one_per_species.tsv   one genome per species (NCBI CheckM species taxid,
                                else the organism taxid): fewest contigs, then earliest
  genomes/<accession>.fna.gz
  seqs.tsv                      accession, sequence id, header (input to topology.py)
"""
import gzip, io, json, os, sys, time, urllib.parse, urllib.request, zipfile

W = sys.argv[1]; SINCE = sys.argv[2] if len(sys.argv) > 2 else '09/05/2025'
API = 'https://api.ncbi.nlm.nih.gov/datasets/v2'
os.makedirs(f'{W}/genomes', exist_ok=True)

def get(url, data=None, timeout=120):
    for t in range(5):
        try:
            req = urllib.request.Request(url, data=data, headers={'Content-Type': 'application/json'} if data else {})
            return urllib.request.urlopen(req, timeout=timeout).read()
        except Exception as e:
            print('retry', e, file=sys.stderr); time.sleep(10 * (t + 1))
    sys.exit(f'failed: {url}')

rep = f'{W}/reports.jsonl'
if not os.path.exists(rep):
    with open(rep + '.tmp', 'w') as out:
        for taxon in (2, 2157):
            tok = None
            while True:
                q = {'filters.assembly_level': 'complete_genome', 'filters.first_release_date': SINCE,
                     'filters.assembly_source': 'genbank', 'page_size': 1000}
                if tok: q['page_token'] = tok
                d = json.loads(get(f'{API}/genome/taxon/{taxon}/dataset_report?' + urllib.parse.urlencode(q)))
                for r in d.get('reports', []): out.write(json.dumps(r) + '\n')
                tok = d.get('next_page_token')
                if not tok: break
                time.sleep(0.4)
    os.rename(rep + '.tmp', rep)
R = [json.loads(l) for l in open(rep)]

best = {}
for r in R:
    ci = r.get('checkm_info', {})
    sp = ci.get('checkm_species_tax_id') or r['organism']['tax_id']
    k = (r['assembly_stats'].get('number_of_contigs', 99), r['assembly_info'].get('release_date', ''))
    if sp not in best or k < best[sp][0]: best[sp] = (k, r)
sel = [v[1] for v in best.values()]
with open(f'{W}/post232_one_per_species.tsv', 'w') as f:
    f.write('accession\tspecies_taxid\torganism\tcontigs\tlength\trelease\tncbi_checkm_compl\tncbi_checkm_cont\n')
    for r in sel:
        ci = r.get('checkm_info', {}); s = r['assembly_stats']
        f.write(f"{r['accession']}\t{ci.get('checkm_species_tax_id', '')}\t{r['organism']['organism_name']}\t"
                f"{s.get('number_of_contigs')}\t{s.get('total_sequence_length')}\t"
                f"{r['assembly_info'].get('release_date')}\t{ci.get('completeness', '')}\t{ci.get('contamination', '')}\n")
print(f'{len(R)} reports, {len(sel)} species')

todo = [r['accession'] for r in sel if not os.path.exists(f"{W}/genomes/{r['accession']}.fna.gz")]
for i in range(0, len(todo), 200):
    body = json.dumps({'accessions': todo[i:i + 200], 'include_annotation_type': ['GENOME_FASTA']}).encode()
    z = zipfile.ZipFile(io.BytesIO(get(f'{API}/genome/download', body, timeout=1800)))
    for n in z.namelist():
        if n.endswith('_genomic.fna'):
            a = '_'.join(os.path.basename(n).split('_')[:2]); p = f'{W}/genomes/{a}.fna.gz'
            with gzip.open(p + '.tmp', 'wb', compresslevel=4) as o: o.write(z.read(n))
            os.rename(p + '.tmp', p)
    print(i + len(todo[i:i + 200]), 'downloaded', flush=True)

with open(f'{W}/seqs.tsv', 'w') as o:
    for r in sel:
        a = r['accession']
        for l in gzip.open(f'{W}/genomes/{a}.fna.gz', 'rt'):
            if l.startswith('>'): h = l[1:].rstrip('\n'); o.write(f'{a}\t{h.split()[0]}\t{h}\n')
