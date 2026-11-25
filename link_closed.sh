#!/usr/bin/env bash
# Link the closed genomes (largest sequence circular) into a run.sh input dir.
#
#   link_closed.sh <genome_topology.tsv> <dir of genomes | gtdb reps database dir> <out dir> [prefix]
#
# A GTDB reps database dir is laid out GCA/123/456/789/<acc>_genomic.fna.gz;
# any other dir holds <acc>.fna.gz. GTDB-Tk refuses genome ids equal to its own
# reference ids, so pass prefix q_ for GTDB representatives.
set -euo pipefail
T=$1 G=$(readlink -f "$2") OUT=$3 P=${4:-}
mkdir -p "$OUT"
awk -F'\t' 'NR>1 && $6=="circular" {print $1}' "$T" | while read -r a; do
    p=$G/$a.fna.gz
    [ -f "$p" ] || p=$G/${a:0:3}/${a:4:3}/${a:7:3}/${a:10:3}/${a}_genomic.fna.gz
    ln -sf "$p" "$OUT/$P$a.fna.gz"
done
echo "$(ls "$OUT" | wc -l) genomes linked in $OUT"
