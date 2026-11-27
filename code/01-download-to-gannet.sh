#!/bin/bash
# Stage the raw reads on gannet. NOT a Hyak job: run it on gannet, or on a
# machine with the gannet volume mounted, inside tmux/screen:
#   bash code/01-download-to-gannet.sh /path/to/gannet/project-ssal-ipn/raw
# Downloads the 36 FASTQ files from Macrogen (85 GB) and checks each against
# the provider's md5sum. Safe to re-run: verified files are skipped and
# partial downloads resume. The links expire two weeks after the report
# date (about 2026-12-10). The directory must be the raw/ folder under
# GANNET_WEB (code/00-config.sh).
set -uo pipefail
dest=${1:?usage: $0 <staging directory>}
links=$(cd "$(dirname "$0")/.." && pwd)/20261126_HN00292532_WGB_Report/assets/spgs/HN00292532_18samples_md5sum_DownloadLink.txt
md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else command md5 -q "$1"; fi; }
mkdir -p "$dest"; cd "$dest" || exit 1
fail=0; n=0
while IFS=$'\t' read -r file size sum url; do
  [ "$file" = File ] && continue
  if grep -qs "^$sum  $file\$" md5-verified.txt; then n=$((n+1)); continue; fi
  wget --continue --tries=5 --no-verbose -O "$file" "$url"
  if [ "$(md5 "$file")" = "$sum" ]; then
    echo "$sum  $file" >> md5-verified.txt; n=$((n+1)); echo "OK   $file"
  else
    echo "FAIL $file (md5 mismatch or incomplete)"; fail=1
  fi
done < <(tr -d '\r' < "$links")
echo "verified $n of 36 files in $dest"
exit $fail
