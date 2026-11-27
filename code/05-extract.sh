#!/bin/bash
#SBATCH --job-name=05-extract
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=2-00:00:00
#SBATCH --array=1-18%6
#SBATCH --output=logs/%x_%A_%a.out
# Methylation extraction from the deduplicated BAM of one library:
# per-context calls, bedGraph/coverage, genome-wide CpG report, and a
# strand-merged CpG coverage file (input for DSS / methylKit / dmrseq).
# Bases ignored at the 5' end of each read default to 3 for R1 and 9 for R2,
# chosen from the M-bias table of the pilot library G1_HK1: R1 is unstable
# over bases 1-3, and leftover adaptase tail inflates non-CpG methylation
# over the first 9 bases of R2 (3.2% falling to the 1.1% baseline). The 3'
# ends are flat. Override with e.g.
#   sbatch --export=ALL,IGNORE_R1=2,IGNORE_R2=12 code/05-extract.sh
#
# When the library is done its outputs go to gannet: the merged CpG coverage
# file is copied (it stays on Hyak for differential methylation), and the
# deduplicated BAM, CpG report, bedGraph, coverage and per-read context
# files are moved. If the push fails they stay on scrubbed until
# code/07-archive-to-gannet.sh is run.
#
# The extractor cannot resume, so a preempted task starts that step again;
# it should fit inside ckpt's 4-5 hour window at this depth. If it does not,
# resubmit on a lab allocation (see code/00-config.sh), for example
#   sbatch --account=coenv --partition=compute --array=3 code/05-extract.sh
set -euo pipefail
source code/00-config.sh
s=$(sample "$SLURM_ARRAY_TASK_ID"); base=$s.deduplicated
mkdir -p "$METH"
if [ ! -f "$METH/$s.done" ]; then
  if [ ! -f "$METH/$s.extract.done" ]; then
    run bismark_methylation_extractor --paired-end --no_overlap --comprehensive \
      --ignore "${IGNORE_R1:-3}" --ignore_r2 "${IGNORE_R2:-9}" \
      --ignore_3prime "${IGNORE_3P_R1:-0}" --ignore_3prime_r2 "${IGNORE_3P_R2:-0}" \
      --gzip --bedGraph --buffer_size 20G \
      --cytosine_report --genome_folder "$GENOME_DIR" \
      --parallel $((SLURM_CPUS_PER_TASK / 3)) -o "$METH" "$DEDUP/$base.bam"
    touch "$METH/$s.extract.done"
  fi
  run coverage2cytosine --merge_CpG --gzip --genome_folder "$GENOME_DIR" \
    --dir "$METH" -o "$s.CpG" "$METH/$base.bismark.cov.gz"

  # CpG coverage summary (strands merged): mean depth over all CpGs in the
  # genome, share of CpGs with at least 1, 3, 5 and 10 reads, and a depth
  # histogram. This is the number that decides what analysis is possible.
  total=$(( $(zcat "$METH/$s.CpG.CpG_report.txt.gz" | wc -l) / 2 ))
  zcat "$METH/$s.CpG.CpG_report.merged_CpG_evidence.cov.gz" |
    awk -v tot="$total" -v s="$s" -v hist="$METH/$s.coverage-hist.tsv" 'BEGIN{OFS="\t"}
      {d=$5+$6; n++; sum+=d; if(d>=3)c3++; if(d>=5)c5++; if(d>=10)c10++; h[d>30?30:d]++}
      END{print "sample","total_cpgs","mean_depth","pct_ge1","pct_ge3","pct_ge5","pct_ge10"
          printf "%s\t%d\t%.2f\t%.1f\t%.1f\t%.1f\t%.1f\n", s, tot, sum/tot, 100*n/tot, 100*c3/tot, 100*c5/tot, 100*c10/tot
          print "depth","cpgs" > hist; print 0, tot-n > hist
          for(i=1;i<=30;i++) print (i==30?"30+":i), h[i]+0 > hist}' > "$METH/$s.coverage.tsv"
  cat "$METH/$s.coverage.tsv"
  touch "$METH/$s.done"
fi
archive_sample "$s" || echo "push to gannet failed; $s stays on scrubbed (run code/07-archive-to-gannet.sh)"
