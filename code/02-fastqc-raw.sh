#!/bin/bash
#SBATCH --job-name=02-fastqc-raw
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=1-00:00:00
#SBATCH --array=1-18
#SBATCH --output=logs/%x_%A_%a.out
# Raw read QC for one library: pull it from gannet, run FastQC, record the
# flowcell and lane, delete the reads. MultiQC across libraries is built by
# code/06-collect-metrics.sh.
set -euo pipefail
source code/00-config.sh
s=$(sample "$SLURM_ARRAY_TASK_ID"); out=$QC/fastqc-raw; mkdir -p "$out"
[ -f "$out/$s.done" ] && { echo "$s already done"; exit 0; }
fetch_raw "$s"
run fastqc -t 2 -o "$out" "$RAW/${s}_1.fastq.gz" "$RAW/${s}_2.fastq.gz"
printf '%s\t%s\n' "$s" "$( (zcat "$RAW/${s}_1.fastq.gz" || true) | head -1 | cut -d: -f1-4)" > "$out/$s.flowcell.tsv"
touch "$out/$s.done"
drop_raw "$s"
to_gannet "qc/fastqc-raw/${s}_*" || echo "push to gannet failed; files stay on scrubbed"
