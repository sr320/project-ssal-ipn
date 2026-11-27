#!/bin/bash
#SBATCH --job-name=02b-multiqc
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=1:00:00
#SBATCH --output=logs/%x_%j.out
# MultiQC across libraries, for whichever FastQC sets exist (raw, trimmed),
# pushed to gannet so the dashboard can link to them. A job rather than a
# login-node command because the container does not start on login nodes.
# Re-run it whenever more libraries have finished.
set -uo pipefail
source code/00-config.sh
if ls "$QC"/fastqc-raw/*_fastqc.zip >/dev/null 2>&1; then
  run multiqc -q -f -o "$QC/multiqc-raw" "$QC/fastqc-raw"
  cat "$QC"/fastqc-raw/*.flowcell.tsv > "$QC/flowcell-lane.tsv"
fi
if ls "$QC"/fastqc-trimmed/*_fastqc.zip >/dev/null 2>&1; then
  run multiqc -q -f -o "$QC/multiqc-trimmed" "$QC/fastqc-trimmed" "$TRIM"
fi
to_gannet qc/multiqc-raw qc/multiqc-trimmed qc/flowcell-lane.tsv || echo "push to gannet failed; reports stay on scrubbed"
