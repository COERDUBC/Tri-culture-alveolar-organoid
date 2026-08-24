#!/usr/bin/env bash
# Step 2 of 3 — CIBERSORTx input files (stages 04-05).
#
#   bash bash_run/step2_cibersortx_prep.sh
#
# Continues the run step 1 created. Stage 05 reads the GSE136831 sparse matrix
# and needs >= 64 GB of RAM. Skip this step on a re-run once the reference TSV
# and the CIBERSORTx result are already named in config.yml.

set -euo pipefail
cd "$(dirname "$0")/.."

export LANG="${LANG:-C.UTF-8}"
export LC_ALL="${LC_ALL:-C.UTF-8}"

Rscript r_scripts/04_cibersortx_prep.R
Rscript r_scripts/05_GSE136831_reference.R

run="$(cat runs/CURRENT_RUN)"
echo
echo "Step 2 complete: runs/${run}"
echo "Upload files: runs/${run}/cibersortx_upload"
echo "Next: run the CIBERSORTx job (see README.md), set"
echo "deconvolution.cibersortx.GSE136831_lineage_result in config_local.yml, then:"
echo "bash bash_run/step3_downstream.sh"
