#!/usr/bin/env bash
# Step 1 of 3 — differential expression and pathway analysis (stages 01-03).
#
#   bash bash_run/step1_de.sh
#
# Creates a new runs/run_YYYY-MM-DD_HH-MM-SS/ folder and records it in
# runs/CURRENT_RUN; steps 2 and 3 continue that same folder.

set -euo pipefail
cd "$(dirname "$0")/.."

# svglite writes UTF-8 one byte at a time under LC_CTYPE=C, which turns the
# figure labels into "TGF-I2" mojibake instead of "TGF-β1".
export LANG="${LANG:-C.UTF-8}"
export LC_ALL="${LC_ALL:-C.UTF-8}"

Rscript r_scripts/01_import.R
Rscript r_scripts/02_de.R
Rscript r_scripts/03_gsea.R

echo
echo "Step 1 complete: runs/$(cat runs/CURRENT_RUN)"
echo "Check that stage 01 printed 'Model: NNNNN genes x 18 samples' above."
echo "Next: bash bash_run/step2_cibersortx_prep.sh"
