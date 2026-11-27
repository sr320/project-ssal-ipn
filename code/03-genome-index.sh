#!/bin/bash
#SBATCH --job-name=03-genome-index
#SBATCH --account=srlab
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
# Download Atlantic salmon Ssal_v3.1 (NCBI RefSeq GCF_905237065.1) and build
# the Bismark bowtie2 bisulfite index. Independent of the reads.
#
# This is the one job sent to the lab partition instead of ckpt: the index
# build cannot resume, and ckpt stops and requeues jobs every 4-5 hours, so
# it could restart from zero indefinitely. To try ckpt anyway:
#   sbatch --partition=ckpt-all code/03-genome-index.sh
# The 30-day scrub applies to the index too; a copy goes to gannet
# (bash code/07-archive-to-gannet.sh genome  if the push here fails).
set -euo pipefail
source code/00-config.sh
mkdir -p "$GENOME_DIR"
[ -f "$GENOME_DIR/index.done" ] && { echo "index already built"; exit 0; }
[ -s "$GENOME_DIR/Ssal_v3.1.fa.gz" ] || wget --no-verbose -O "$GENOME_DIR/Ssal_v3.1.fa.gz" "$GENOME_URL"
run bismark_genome_preparation --verbose --parallel $((SLURM_CPUS_PER_TASK / 2)) "$GENOME_DIR"
touch "$GENOME_DIR/index.done"
[ "$GENOME_DIR" = "$WORK/genome" ] && { to_gannet genome || echo "push to gannet failed; index stays on scrubbed only"; }
