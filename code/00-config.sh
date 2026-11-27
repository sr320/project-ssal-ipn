# Shared settings for every Hyak (klone) job. Sourced, not run.
# Submit all jobs from the repository root:  sbatch code/NN-name.sh
#
# Where jobs run. Scripts default to the checkpoint partition. If ckpt is
# not enough (a library keeps getting preempted, or the queue is slow),
# resubmit with one of the lab's own allocations, which are not preempted:
#   sbatch --account=srlab --partition=cpu-g2-mem2x ...   (all of it)
#   sbatch --account=coenv --partition=compute ...
#   sbatch --account=coenv --partition=cpu-g2 ...          (up to 28 CPUs)
# Command-line options override the #SBATCH lines in the script.
#
# Layout: gannet holds the data, Hyak only computes.
#   - Raw reads are staged on gannet. Each job pulls the one library it
#     needs over https, works on scrubbed, and deletes the reads afterwards.
#   - Finished files are pushed back to gannet with rsync, and files that
#     are no longer needed on Hyak are removed there once the push succeeds.
#   - The dashboard links to the files on gannet.

# --- gannet project directory: /volume1/v1_web/project-ssal-ipn ------------
# It will hold  raw/  (the 36 FASTQ files) and  output/  (results).
# Web address of that directory (read by Hyak jobs and the dashboard):
GANNET_WEB=${GANNET_WEB:-https://gannet.fish.washington.edu/v1_web/project-ssal-ipn}
# The same directory as an rsync/ssh destination:
GANNET_SSH=${GANNET_SSH:-$USER@gannet.fish.washington.edu:/volume1/v1_web/project-ssal-ipn}
# ---------------------------------------------------------------------------

REPO=${SLURM_SUBMIT_DIR:-$(pwd)}

# Scratch space: 5 TB per user, files deleted 30 days after creation.
WORK=${WORK:-/gscratch/scrubbed/$USER/project-ssal-ipn}
RAW=$WORK/raw                      # transient, one library per job
TRIM=$WORK/trim                    # reads transient; logs kept
ALN=$WORK/bismark                  # chunk BAMs transient; reports kept
DEDUP=$WORK/dedup                  # deduplicated BAMs, until extraction
METH=$WORK/methylation             # per-library methylation calls
QC=$WORK/qc
GENOME_DIR=${GENOME_DIR:-$WORK/genome}   # Ssal_v3.1.fa.gz + Bisulfite_Genome/

SAMPLES=$REPO/data/samples.tsv
LINKS=$REPO/20261126_HN00292532_WGB_Report/assets/spgs/HN00292532_18samples_md5sum_DownloadLink.txt
GENOME_URL=https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/905/237/065/GCF_905237065.1_Ssal_v3.1/GCF_905237065.1_Ssal_v3.1_genomic.fna.gz

# Lab bioinformatics container (newest build unless SIF is set).
SIF=${SIF:-$(ls -t /gscratch/srlab/containers/srlab-bioinformatics-container*.sif 2>/dev/null | head -1)}
run() {  # run a program inside the container
  apptainer exec --bind /mmfs1/home/ --bind /mmfs1/gscratch/ --bind /gscratch/ "$SIF" "$@"
}

# Sample name for an array task (1-18), or all sample names.
sample() { awk -v n="$1" 'NR==n+1{print $1}' "$SAMPLES"; }
all_samples() { awk 'NR>1{print $1}' "$SAMPLES"; }

# Pull both FASTQ files of one library from gannet into $RAW and check them
# against Macrogen's md5sums. Resumes partial files; skips verified ones.
# RAW_SOURCE=macrogen fetches from the Macrogen links instead (for a pilot
# before gannet is set up; the links expire about 2026-12-10).
fetch_raw() {
  local s=$1 r f want url
  mkdir -p "$RAW"
  for r in 1 2; do
    f=${s}_$r.fastq.gz
    [ -f "$RAW/$f.ok" ] && [ -s "$RAW/$f" ] && continue
    want=$(tr -d '\r' < "$LINKS" | awk -F'\t' -v f="$f" '$1==f{print $3}')
    url=$GANNET_WEB/raw/$f
    [ "${RAW_SOURCE:-gannet}" = macrogen ] && url=$(tr -d '\r' < "$LINKS" | awk -F'\t' -v f="$f" '$1==f{print $4}')
    wget --continue --tries=5 --no-verbose -O "$RAW/$f" "$url" || return 1
    if [ "$(md5sum "$RAW/$f" | cut -d' ' -f1)" = "$want" ]; then
      touch "$RAW/$f.ok"
    else
      echo "md5 mismatch for $f; removing it" >&2; rm -f "${RAW:?}/$f"; return 1
    fi
  done
}
drop_raw() { rm -f "${RAW:?}/${1:?}"_[12].fastq.gz "${RAW:?}/${1:?}"_[12].fastq.gz.ok; }

# Push files to gannet, keeping their path relative to $WORK under output/.
#   to_gannet [--move] <path or glob relative to $WORK> ...
# --move deletes each file on Hyak after it has transferred. Patterns that
# match nothing are skipped. Inside jobs this runs without a password prompt
# (GANNET_BATCH=1), so it needs an ssh key for gannet; if the push fails the
# files stay on scrubbed and code/07-archive-to-gannet.sh picks them up.
GANNET_BATCH=${GANNET_BATCH:-1}
to_gannet() {
  local opts=(--archive --relative) files=() p q
  [ "${1:-}" = --move ] && { opts+=(--remove-source-files); shift; }
  [ "$GANNET_BATCH" = 1 ] && opts+=(-e "ssh -o BatchMode=yes -o ConnectTimeout=30")
  cd "$WORK" || return 1
  for p in "$@"; do for q in $p; do [ -e "$q" ] && files+=("./$q"); done; done   # $p unquoted: expand globs here
  [ ${#files[@]} -eq 0 ] && { cd "$OLDPWD"; return 0; }
  rsync "${opts[@]}" "${files[@]}" "$GANNET_SSH/output/"
  local rc=$?; cd "$OLDPWD"; return $rc
}

# Archive one library once its methylation extraction is done. The merged
# CpG coverage file stays on Hyak as well (differential methylation reads
# it); the deduplicated BAM and the other large outputs move to gannet.
archive_sample() {
  local s=$1
  [ -f "$METH/$s.done" ] || { echo "$s: extraction not finished, nothing archived"; return 0; }
  [ -f "$METH/$s.archived" ] && return 0
  to_gannet methylation/$s.CpG.*merged_CpG_evidence.cov.gz methylation/$s.coverage*.tsv \
            methylation/$s.deduplicated*report.txt methylation/$s.deduplicated.M-bias.txt \
            methylation/mbias/$s.* dedup/$s.deduplication_report.txt \
            bismark/${s}_1.chunk*_PE_report.txt trim/$s.pass?.log &&
  to_gannet --move dedup/$s.deduplicated.bam \
            methylation/$s.deduplicated.CpG_report.txt.gz methylation/$s.CpG.CpG_report.txt.gz \
            methylation/$s.deduplicated.bismark.cov.gz methylation/$s.deduplicated.bedGraph.gz \
            methylation/*_context_$s.deduplicated.txt.gz &&
  touch "$METH/$s.archived" && echo "$s archived to gannet"
}

module load apptainer 2>/dev/null || true
mkdir -p "$WORK"
