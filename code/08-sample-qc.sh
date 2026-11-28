#!/bin/bash
#SBATCH --job-name=08-sample-qc
#SBATCH --account=srlab
#SBATCH --partition=ckpt-all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=logs/%x_%j.out
# Sample-level QC across all libraries, on windows rather than single CpGs
# (libraries have about 2x per CpG, so single-site values are mostly noise).
#
#   1. Per library, sum methylated and total read counts of the strand-merged
#      CpG coverage file over fixed windows (TILE bp, default 1000).
#   2. code/08-sample-qc.R: how many windows pass each coverage filter,
#      per-library and per-chromosome methylation, correlation, clustering
#      and PCA across libraries, with tables and plots.
#   3. Push the tables and plots to gannet (output/qc/sample-qc/).
#
# Coverage files missing from scrubbed are fetched from gannet. Each step
# leaves a marker, so a requeued job continues where it stopped.
#   sbatch code/08-sample-qc.sh
#   sbatch --export=ALL,TILE=500 code/08-sample-qc.sh
set -uo pipefail
source code/00-config.sh
TILE=${TILE:-1000}
OUT=$QC/sample-qc
TILES=$OUT/tiles-$TILE
mkdir -p "$TILES" logs

tile_one() {
  local s=$1 cov=$METH/$1.CpG.CpG_report.merged_CpG_evidence.cov.gz
  [ -f "$TILES/$s.done" ] && return 0
  if [ ! -s "$cov" ]; then
    wget -q -O "$cov.part" "$GANNET_WEB/output/methylation/$(basename "$cov")" \
      && mv "$cov.part" "$cov" || { echo "$s: no coverage file on scrubbed or gannet"; return 1; }
  fi
  # columns in: chr start end pct meth unmeth; out: chr window_start meth total cpgs
  zcat "$cov" | awk -v w="$TILE" 'BEGIN { OFS = "\t" }
    { k = $1 OFS int(($2 - 1) / w) * w; m[k] += $5; t[k] += $5 + $6; n[k]++ }
    END { for (k in t) print k, m[k], t[k], n[k] }' | gzip > "$TILES/$s.tsv.gz.part" \
    && mv "$TILES/$s.tsv.gz.part" "$TILES/$s.tsv.gz" && touch "$TILES/$s.done" \
    && echo "$s: windows written"
}
export -f tile_one
export METH TILES TILE GANNET_WEB

all_samples | xargs -P "${SLURM_CPUS_PER_TASK:-4}" -I{} bash -c 'tile_one {}'
missing=0
for s in $(all_samples); do [ -f "$TILES/$s.done" ] || { echo "missing windows for $s"; missing=1; }; done
[ "$missing" = 0 ] || exit 1

run Rscript code/08-sample-qc.R "$SAMPLES" "$TILES" "$OUT" "$TILE" || exit 1
touch "$OUT/sample-qc.done"
to_gannet "qc/sample-qc/*.tsv" "qc/sample-qc/*.pdf" "qc/sample-qc/*.png" \
  || echo "push to gannet failed; results stay in $OUT"
