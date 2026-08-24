#!/usr/bin/env bash
# Step 3 of 3 — deconvolution and cross-talk (stages 06-08).
#
#   bash bash_run/step3_downstream.sh
#
# Continues the run step 1 created. Set
# deconvolution.cibersortx.GSE136831_lineage_result in config_local.yml first.

set -euo pipefail
cd "$(dirname "$0")/.."

export LANG="${LANG:-C.UTF-8}"
export LC_ALL="${LC_ALL:-C.UTF-8}"

Rscript r_scripts/06_deconvolution.R
Rscript r_scripts/07_lr.R
# Stage 08 prints a skip message when its NicheNet files are blank or TO_FILL.
Rscript r_scripts/08_nichenet.R

echo
echo "Step 3 complete: runs/$(cat runs/CURRENT_RUN)"
