#!/usr/bin/env Rscript
# One-time package installation for the analysis pipeline.

cran <- c("yaml", "dplyr", "tidyr", "readr", "tibble",
          "ggplot2", "ggrepel", "pheatmap", "RColorBrewer", "matrixStats",
          "svglite", "msigdbr", "here", "png", # export alpha-channel tests
          "circlize")   # required for Stage 07/08 L-R chord diagrams
new_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_cran)) install.packages(new_cran, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

bioc <- c("tximport", "limma", "edgeR",
          "AnnotationDbi", "fgsea",
          "org.Hs.eg.db")  # human; swap to org.Mm.eg.db etc. for other species (match config.yml)
new_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(new_bioc)) BiocManager::install(new_bioc, ask = FALSE, update = FALSE)

# Note: CIBERSORTx runs through its website; no local CIBERSORT package is
# required.
# Stage 07 cross-cell L-R uses a curated CSV (connectomeDB2025, set crosstalk$lr_csv
# in config.yml) — no L-R package or live database dependency.
# Optional: NicheNet ligand-activity prioritization (Stage 08). Install it from
# GitHub only when it is absent, so rerunning setup does not reinstall it.
if (!requireNamespace("nichenetr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("saeyslab/nichenetr")
}

message("Setup complete.")
