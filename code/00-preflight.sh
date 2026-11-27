#!/bin/bash
# Check the container, the gannet addresses and storage before submitting
# jobs. Run on a klone login node:  bash code/00-preflight.sh
source code/00-config.sh
echo "container: ${SIF:-NOT FOUND}"
echo "work dir:  $WORK"
echo "core dir:  $CORE  (genome index: $([ -f "$GENOME_DIR/index.done" ] && echo built || echo not built yet))"
ok=1
for p in fastqc multiqc bismark bismark_genome_preparation \
         deduplicate_bismark bismark_methylation_extractor coverage2cytosine \
         bowtie2 samtools; do
  if run bash -c "command -v $p" >/dev/null 2>&1; then
    printf '  ok       %s\n' "$p"
  else
    printf '  MISSING  %s\n' "$p"; ok=0
  fi
done
if v=$(cutadapt_run --version 2>/dev/null); then
  printf '  ok       cutadapt %s\n' "$v"
else
  echo "  MISSING  cutadapt  (install it with  bash code/00-install-cutadapt.sh)"; ok=0
fi
first=$(all_samples | head -1)_1.fastq.gz
if wget --spider --quiet "$GANNET_WEB/raw/$first"; then
  echo "  ok       raw reads on gannet ($GANNET_WEB/raw)"
else
  echo "  MISSING  $GANNET_WEB/raw/$first"
  echo "           stage the reads with code/01-download-to-gannet.sh, or fix GANNET_WEB in code/00-config.sh"; ok=0
fi
host=${GANNET_SSH%%:*}
if ssh -o BatchMode=yes -o ConnectTimeout=15 "$host" true 2>/dev/null; then
  echo "  ok       passwordless ssh to $host (jobs can push results themselves)"
else
  echo "  note     no passwordless ssh to $host: jobs will leave results on scrubbed;"
  echo "           push them with  bash code/07-archive-to-gannet.sh  (asks for a password), or add an ssh key"
fi
[ "$ok" = 1 ] && echo "ready" || echo "not ready: fix the items marked MISSING"
