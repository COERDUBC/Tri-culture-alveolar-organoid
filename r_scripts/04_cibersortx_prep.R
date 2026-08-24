#!/usr/bin/env Rscript
# 04 — Write the TSV files to upload to the CIBERSORTx web tool.
# CIBERSORTx is not run here: it is a hosted tool, so this stage prepares inputs
# and Stage 06 reads the results you download back.
#
# Upload format: tab-separated, first column the gene symbol, one column per
# sample / cell type, no quoting.

source("r_scripts/common.R")

m <- readRDS(obj_path("model.rds"))
up <- file.path(OUT, "cibersortx_upload")
dir.create(up, recursive = TRUE, showWarnings = FALSE)

write_ref <- function(mat, file) {
  df <- data.frame(GeneSymbol = rownames(mat), mat, check.names = FALSE)
  write.table(df, file.path(up, file), sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = TRUE)
  message("  ", file, ": ", nrow(mat), " genes x ", ncol(mat), " columns")
}

# Mixture: the tri-culture bulk, CPM, symbol-collapsed.
mixture <- edgeR::cpm(to_symbol(m$dge$counts), log = FALSE)
write_ref(round(mixture, 4), "bulk_mixture_cpm.tsv")

cat(sprintf('
Next steps (https://cibersortx.stanford.edu):
  1. Upload from %s
       Mixture       : bulk_mixture_cpm.tsv
       Reference     : GSE136831_reference_donor_cap<k>.tsv (from Stage 05)
  2. Use the reference to create a five-state signature matrix, then run
     "Impute Cell Fractions" with bulk_mixture_cpm.tsv and
       permutations          = %s
       quantile normalisation = %s   (off is correct for RNA-seq)
  3. Save the downloaded GSE136831 result CSV, then point config.yml at it:
       deconvolution:
         cibersortx:
           GSE136831_lineage_result: "<path to downloaded CSV>"
  4. Rscript r_scripts/06_deconvolution.R
',
  up,
  cfg$deconvolution$permutations %||% 1000,
  if (isTRUE(cfg$deconvolution$quantile_normalization)) "on" else "off"))

message("04 done -> ", up)
