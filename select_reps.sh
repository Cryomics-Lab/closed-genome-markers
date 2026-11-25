#!/usr/bin/env bash
# GTDB species representatives assembled to NCBI "Complete Genome", with the
# FASTA headers of every sequence in them.
#
#   select_reps.sh <dir with bac120/ar53_metadata_rNNN.tsv.gz> <unpacked gtdb_genomes_reps_rNNN/database> <workdir>
#
# Writes <workdir>/complete_reps.tsv (accession, contig count) and
# <workdir>/seqs.tsv (accession, sequence id, full header), the input to topology.py.
set -euo pipefail
META=$(readlink -f "$1") G=$(readlink -f "$2") W=$3
mkdir -p "$W/hdr"; cd "$W"
zcat "$META"/bac120_metadata_r*.tsv.gz "$META"/ar53_metadata_r*.tsv.gz | awk -F'\t' '
    NR==1 {for (i=1; i<=NF; i++) c[$i]=i; next}
    $1=="accession" {next}
    $c["gtdb_representative"]=="t" && $c["ncbi_assembly_level"]=="Complete Genome" {
        print substr($1,4) "\t" $c["contig_count"]}' > complete_reps.tsv
echo "$(wc -l < complete_reps.tsv) complete representatives"
# One header file per genome: parallel writers must not share an output file.
cut -f1 complete_reps.tsv | G=$G xargs -P 16 -I{} bash -c '
    a={}; p=$G/${a:0:3}/${a:4:3}/${a:7:3}/${a:10:3}/${a}_genomic.fna.gz
    [ -f "$p" ] || { echo "$a missing" >&2; exit 0; }
    zcat "$p" | grep "^>" | sed "s/^>//" | awk -v a="$a" "{print a\"\t\"\$1\"\t\"\$0}" > hdr/$a.tsv'
cat hdr/*.tsv > seqs.tsv
echo "$(wc -l < seqs.tsv) sequences"
