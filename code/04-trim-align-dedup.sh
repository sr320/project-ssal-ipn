#!/bin/bash
#SBATCH --job-name=04-trim-align
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --time=3-00:00:00
#SBATCH --array=1-18%6
#SBATCH --output=logs/%x_%A_%a.out
# One library, end to end, streaming from gannet:
#   fetch raw reads -> trim -> FastQC -> split into chunks -> Bismark per
#   chunk -> deduplicate all chunks together -> M-bias -> delete raw reads,
#   trimmed reads and chunk BAMs.
# Only the deduplicated BAM (a few GB) and the reports stay on scrubbed;
# code/05-extract.sh moves the BAM to gannet when it has been used.
#
# Built for ckpt, which stops and requeues jobs every 4-5 hours: each stage
# and each alignment chunk leaves a marker, so a requeued task picks up
# where it stopped and repeats at most one chunk. If a library keeps
# getting preempted, send it to a lab allocation (see code/00-config.sh):
#   sbatch --partition=cpu-g2-mem2x --array=7 code/04-trim-align-dedup.sh
#   sbatch --account=coenv --partition=compute --array=7 code/04-trim-align-dedup.sh
#   sbatch --account=coenv --partition=cpu-g2 --cpus-per-task=28 --array=7 code/04-trim-align-dedup.sh
#
# Settings (override with  sbatch --export=ALL,NAME=value ...):
#   CLIP    bp clipped from the R2 5' end and R1 3' end (adaptase tail).
#           15 cleared it in the earlier salmon WGBS project with this kit;
#           confirm against the raw FastQC per-base content plots.
#   MAXINS  maximum fragment length for a valid pair. Bismark's default of
#           500 discarded 10-17% of pairs in that project; 1000 kept them.
#   SUBSET  pilot mode: align only the first SUBSET pairs, print the mapping
#           efficiency, and stop. Use it to compare MAXINS values:
#             sbatch --array=1 --export=ALL,SUBSET=200000,MAXINS=500  code/04-trim-align-dedup.sh
#             sbatch --array=1 --export=ALL,SUBSET=200000,MAXINS=1000 code/04-trim-align-dedup.sh
#   CHUNK   read pairs per alignment chunk (default 10 million).
#   PARALLEL, P   Bismark instances and bowtie2 threads each (about
#           PARALLEL*2*P cores and 15-20 GB RAM per instance).
#   KEEP_TRIMMED=1   move the trimmed reads to gannet instead of deleting
#           them (they can be regenerated from the raw reads).
set -euo pipefail
source code/00-config.sh
s=$(sample "$SLURM_ARRAY_TASK_ID"); t=$SLURM_CPUS_PER_TASK
CLIP=${CLIP:-15}; MAXINS=${MAXINS:-1000}; CHUNK=${CHUNK:-10000000}
PARALLEL=${PARALLEL:-7}; P=${P:-2}; SUBSET=${SUBSET:-}
mkdir -p "$TRIM" "$ALN" "$DEDUP" "$METH/mbias" "$QC/fastqc-trimmed"
[ -f "$DEDUP/$s.done" ] && [ -z "$SUBSET" ] && { echo "$s already finished"; exit 0; }
[ -f "$GENOME_DIR/index.done" ] || { echo "genome index not built (code/03-genome-index.sh)"; exit 1; }
r1=$TRIM/${s}_1.trimmed.fastq.gz; r2=$TRIM/${s}_2.trimmed.fastq.gz

# 1. Fetch and trim.
#    pass 1: Illumina adapters, 2-colour poly-G/quality trim (q20), clip R2 5' end
#    pass 2: clip R1 3' end; drop pairs shorter than 20 bp
if [ ! -f "$TRIM/$s.trim.done" ]; then
  fetch_raw "$s"
  cutadapt_run -j "$t" --interleaved \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    --nextseq-trim 20 -U "$CLIP" \
    "$RAW/${s}_1.fastq.gz" "$RAW/${s}_2.fastq.gz" 2> "$TRIM/$s.pass1.log" |
  cutadapt_run -j "$t" --interleaved -u "-$CLIP" -m 20 -o "$r1" -p "$r2" - \
    > "$TRIM/$s.pass2.log" 2>&1
  touch "$TRIM/$s.trim.done"
  drop_raw "$s"
fi

# Pilot: align a subset and stop.
if [ -n "$SUBSET" ]; then
  out=$WORK/pilot/maxins$MAXINS; mkdir -p "$out"
  for r in 1 2; do
    (zcat "$TRIM/${s}_$r.trimmed.fastq.gz" || true) | head -n $((SUBSET * 4)) | gzip -1 > "$out/${s}_$r.pilot.fq.gz"
  done
  run bismark --genome "$GENOME_DIR" -1 "$out/${s}_1.pilot.fq.gz" -2 "$out/${s}_2.pilot.fq.gz" \
    --parallel 4 -p "$P" --maxins "$MAXINS" --temp_dir "$out/tmp-$s" -o "$out"
  grep -E "Mapping efficiency|C methylated in" "$out/${s}_1.pilot_bismark_bt2_PE_report.txt"
  exit 0
fi

# 2. FastQC on trimmed reads.
if [ ! -f "$QC/fastqc-trimmed/$s.done" ]; then
  run fastqc -t 2 -o "$QC/fastqc-trimmed" "$r1" "$r2"
  touch "$QC/fastqc-trimmed/$s.done"
fi

# 3. Split into chunks of CHUNK pairs (R1 and R2 stay in step).
cdir=$TRIM/chunks-$s
if [ ! -f "$cdir/split.done" ]; then
  rm -rf "${cdir:?}"; mkdir -p "$cdir"
  for r in 1 2; do
    zcat "$TRIM/${s}_$r.trimmed.fastq.gz" |
      split -l $((CHUNK * 4)) -d -a 2 --additional-suffix=.fq --filter='gzip -1 > $FILE.gz' - "$cdir/${s}_$r.chunk"
  done
  touch "$cdir/split.done"
fi

# 4. Align each chunk (directional, paired-end).
for c1 in "$cdir/${s}"_1.chunk*.fq.gz; do
  k=$(basename "$c1" .fq.gz); k=${k##*chunk}
  [ -f "$ALN/$s.chunk$k.done" ] && continue
  run bismark --genome "$GENOME_DIR" -1 "$c1" -2 "$cdir/${s}_2.chunk$k.fq.gz" \
    --parallel "$PARALLEL" -p "$P" --maxins "$MAXINS" \
    --temp_dir "$ALN/tmp-$s" -o "$ALN"
  rm -rf "${ALN:?}/tmp-$s"
  touch "$ALN/$s.chunk$k.done"
done

# 5. Remove PCR duplicates across all chunks of the library.
ddir=$DEDUP/tmp-$s; rm -rf "${ddir:?}"; mkdir -p "$ddir"
run deduplicate_bismark --paired --bam --multiple --output_dir "$ddir" \
  "$ALN/${s}"_1.chunk*_bismark_bt2_pe.bam
mv "$ddir"/*deduplicated.bam "$DEDUP/$s.deduplicated.bam"
mv "$ddir"/*deduplication_report.txt "$DEDUP/$s.deduplication_report.txt"
rmdir "$ddir"

# 6. M-bias table (methylation by read position), to choose the ignore
#    settings for code/05-extract.sh.
run bismark_methylation_extractor --mbias_only --paired-end --no_overlap \
  --parallel $((t / 3)) -o "$METH/mbias" "$DEDUP/$s.deduplicated.bam"

# 7. Finished: push the small files to gannet, drop the bulky intermediates.
touch "$DEDUP/$s.done"
rm -rf "${cdir:?}" "${ALN:?}/${s}"_1.chunk*_bismark_bt2_pe.bam
to_gannet "qc/fastqc-trimmed/${s}_*" "trim/$s.pass?.log" "bismark/${s}_1.chunk*_PE_report.txt" \
          "dedup/$s.deduplication_report.txt" "methylation/mbias/$s.*" \
  || echo "push to gannet failed; reports stay on scrubbed"
if [ "${KEEP_TRIMMED:-0}" = 1 ] && to_gannet --move "trim/${s}_?.trimmed.fastq.gz"; then
  echo "trimmed reads moved to gannet"
else
  rm -f "$r1" "$r2"
fi
echo "$s finished $(date)"
