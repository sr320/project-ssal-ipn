#!/bin/bash
# Catch-up push from scrubbed to gannet, for anything the jobs could not
# send themselves (no ssh key on the compute node, gannet unreachable).
# Run on a klone login node inside tmux; it asks for your gannet password
# unless ssh keys are set up.
#   bash code/07-archive-to-gannet.sh
# For each library whose extraction is finished: the merged CpG coverage
# file is copied, and the deduplicated BAM and other large outputs are moved
# (deleted on Hyak after they transfer). Then run code/06-collect-metrics.sh
# so the dashboard links appear.
set -uo pipefail
source code/00-config.sh
GANNET_BATCH=0
to_gannet qc "trim/*.log" "bismark/*_PE_report.txt" "dedup/*deduplication_report.txt" methylation/mbias
for s in $(all_samples); do archive_sample "$s"; done
echo "archived: $(ls "$METH"/*.archived 2>/dev/null | wc -l) of 18 libraries"
