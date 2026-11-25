# Marker genes and CheckM2 on closed genomes

Where do GTDB-Tk's marker genes and CheckM2 fall short on genomes that are
known to be complete? A closed, circular genome should have every
single-copy marker exactly once, 100% completeness and 0% contamination, so
any departure from that is an error of the tool or its marker set, not of the
genome. We score two sets:

- **reps**: the 10,737 GTDB R232 species representatives whose NCBI assembly
  is "Complete Genome" and whose largest sequence is circular (the best case:
  GTDB built its marker sets and CheckM2 trained on genomes like these);
- **post232**: 4,744 complete genomes released to GenBank after R232's NCBI
  cutoff (2025-09-04), one per species, largest sequence circular (held out).

No genome is altered: no fragmentation, no added contamination. Only the
GTDB-Tk steps that use marker genes run (`identify`, `align`); nothing is
classified.

Write-up: https://cryomics.org/reports/closed-genome-markers/

## Pipeline

| Step | Script | Output |
|---|---|---|
| Complete GTDB representatives and their sequence headers | `select_reps.sh` | `complete_reps.tsv`, `seqs.tsv` |
| Post-R232 complete genomes, one per species, downloaded | `select_post232.py` | `post232_one_per_species.tsv`, `genomes/`, `seqs.tsv` |
| NCBI topology of every sequence (closed = largest sequence circular) | `topology.py` | `genome_topology.tsv` |
| Link the closed genomes into an input dir | `link_closed.sh` | `in/` |
| CheckM2 predict; GTDB-Tk identify and align | `run.sh` | `out/` |
| GTDB lineage for post-R232 genomes, from their NCBI names | `ncbi_to_gtdb_lineage.py` | `post232_lineage.tsv` |
| Score each genome against its domain's marker set | `markers.py` | `<set>.genomes.tsv`, `<set>.markers.tsv` |
| Small CSVs for the report's figures | `report_data.py` | `summary.csv`, `dist.csv`, ... |

```bash
R=/path/to/gtdb/release232        # bac120/ar53_metadata_r232.tsv.gz, gtdb_genomes_reps_r232 unpacked
./select_reps.sh $R $R/gtdb_genomes_reps_r232/database reps
python3 topology.py reps
./link_closed.sh reps/genome_topology.tsv $R/gtdb_genomes_reps_r232/database reps/in q_
./run.sh reps/in reps/out

python3 select_post232.py post232
python3 topology.py post232
./link_closed.sh post232/genome_topology.tsv post232/genomes post232/in
./run.sh post232/in post232/out

zcat $R/bac120_taxonomy_r232.tsv.gz $R/ar53_taxonomy_r232.tsv.gz | sed 's/^[A-Z][A-Z]_//' > gtdb_r232_taxonomy.tsv
python3 ncbi_to_gtdb_lineage.py $R/bac120_metadata_r232.tsv.gz $R/ar53_metadata_r232.tsv.gz \
    post232/post232_one_per_species.tsv > post232_lineage.tsv
python3 markers.py reps/out gtdb_r232_taxonomy.tsv reps
python3 markers.py post232/out post232_lineage.tsv post232
python3 report_data.py . reps/out/gtdbtk/align/intermediate_results report_data
```

The `q_` prefix is needed because GTDB-Tk refuses input genomes whose ids
equal its reference genome ids. NCBI requests carry `tool=gtdb_topology` and
stay under 3 per second.

## Versions

GTDB R232 (GTDB-Tk reference data r232), GTDB-Tk 2.8.0, CheckM2 1.1.0 with
its `uniref100.KO.1` DIAMOND database, NCBI Datasets API v2 and E-utilities
(queried 2026-11-24).

## Results

`results/` holds the per-genome and per-marker tables (`<set>.genomes.tsv`,
`<set>.markers.tsv`), each set's `genome_topology.tsv`, the post-R232 picks
and their GTDB lineages.

## Licence

MIT
