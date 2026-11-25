#!/usr/bin/env bash
# CheckM2 and GTDB-Tk marker identification/alignment on clean closed genomes.
#   run.sh <genome dir of *.fna.gz> <out dir>
# Only the GTDB-Tk steps that use marker genes (identify, align) run; no
# classification. Tool paths default to this lab's installs; override with
# CHECKM2, CHECKM2_DB, GTDBTK_BIN and GTDBTK_DATA_PATH. A finished step is skipped.
set -euo pipefail
IN=$(readlink -f "$1"); OUT=$(readlink -f "$2"); mkdir -p "$OUT"
T=${THREADS:-30}
CK=${CHECKM2:-/opt/miniforge3/envs/checkm2/bin/checkm2}
CKDB=${CHECKM2_DB:-/data/scratch/refdbs/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd}
GT=${GTDBTK_BIN:-/home/cryomics/.local/share/mamba/envs/gtdbtk/bin}
export GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-/data/scratch/refdbs/gtdbtk_db}
if [ ! -s "$OUT/checkm2/quality_report.tsv" ]; then
    echo "[$(date +%T)] checkm2"
    $CK predict -i "$IN" -x fna.gz -o "$OUT/checkm2" -t $T --force \
        --database_path "$CKDB" \
        > "$OUT/checkm2.log" 2>&1
fi
if [ ! -d "$OUT/gtdbtk/identify" ]; then
    echo "[$(date +%T)] gtdbtk identify"
    PATH=$GT:$PATH $GT/gtdbtk identify --genome_dir "$IN" -x fna.gz --out_dir "$OUT/gtdbtk" \
        --cpus $T > "$OUT/gtdbtk_identify.log" 2>&1
fi
if [ ! -d "$OUT/gtdbtk/align" ]; then
    echo "[$(date +%T)] gtdbtk align"
    PATH=$GT:$PATH $GT/gtdbtk align --identify_dir "$OUT/gtdbtk" --out_dir "$OUT/gtdbtk" \
        --cpus $T > "$OUT/gtdbtk_align.log" 2>&1
fi
echo "[$(date +%T)] done"
