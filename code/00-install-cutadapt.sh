#!/bin/bash
# The lab container has no cutadapt. Install it once into lab group storage
# ($CORE/tools/cutadapt, a few MB) using the container's own Python, so it
# persists and runs inside the container like everything else.
# Run on a compute node (inside salloc):  bash code/00-install-cutadapt.sh
set -euo pipefail
source code/00-config.sh
dest=$CORE/tools/cutadapt
py=$(run bash -c 'echo "$(dirname "$(command -v pip)")/python"')
echo "installing cutadapt with $py into $dest"
mkdir -p "$dest"
run "$py" -m pip install --no-cache-dir --upgrade --target "$dest" cutadapt
echo "$py" > "$dest/PYTHON"
echo "installed: cutadapt $(cutadapt_run --version)"
