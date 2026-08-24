#!/usr/bin/env Rscript
# 01 — Import Salmon counts, fit the limma-voom model, sample QC (PCA, distances).
# Batch is a covariate in the model; removeBatchEffect is used for PLOTS ONLY.

# Stage 01 is the start of a new analysis. If OUT was explicitly supplied,
# common.R honours it; otherwise this option creates a fresh timestamped run.
options(fibres.start_new_run = TRUE)
source("r_scripts/common.R")
suppressPackageStartupMessages({
  library(tximport); library(edgeR); library(limma); library(matrixStats)
})

# --- samples ---------------------------------------------------------------
samples <- read_csv(cfg$input$metadata, show_col_types = FALSE) |> as.data.frame()
cond <- cfg$design$condition_column
samples <- samples[!is.na(samples[[cond]]) & samples[[cond]] != "NA", ]
rownames(samples) <- samples$sample
configured_conditions <- names(condition_label_map())
observed_conditions <- unique(as.character(samples[[cond]]))
if (!identical(configured_conditions[[1]], cfg$design$reference_level))
  stop("The first display_labels$conditions entry must be the design reference ",
       "level (", cfg$design$reference_level, ").", call. = FALSE)
if (!setequal(observed_conditions, configured_conditions))
  stop("Observed conditions do not match display_labels$conditions. ",
       set_diff_text(configured_conditions, observed_conditions),
       call. = FALSE)
samples[[cond]] <- factor(samples[[cond]], levels = configured_conditions)
for (v in setdiff(all.vars(as.formula(cfg$design$formula)), cond))
  samples[[v]] <- factor(samples[[v]])

# --- counts ----------------------------------------------------------------
files <- file.path(cfg$input$salmon_dir, rownames(samples), "quant.sf")
names(files) <- rownames(samples)
if (any(!file.exists(files))) stop("Missing quant.sf: ",
                                   paste(names(files)[!file.exists(files)], collapse = ", "))

t2g <- read_tsv(cfg$input$tx2gene, col_names = FALSE, show_col_types = FALSE)
t2g[[1]] <- sub("\\.\\d+$", "", t2g[[1]])          # strip transcript version
txi <- tximport(files, type = "salmon", tx2gene = t2g[, 1:2], ignoreTxVersion = TRUE)

# gene_id -> symbol: tx2gene column 3, then OrgDb for whatever is left.
sym <- if (ncol(t2g) >= 3) {
  m <- distinct(t2g[, 2:3]); setNames(m[[2]], m[[1]])
} else setNames(character(0), character(0))
sym[sym == "" | sym == names(sym)] <- NA
gaps <- setdiff(rownames(txi$counts), names(sym)[!is.na(sym)])
if (length(gaps) && requireNamespace(cfg$project$orgdb, quietly = TRUE)) {
  db <- getExportedValue(cfg$project$orgdb, cfg$project$orgdb)
  hit <- suppressMessages(tryCatch(AnnotationDbi::mapIds(
    db, keys = gaps, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"),
    error = function(e) NULL))
  if (!is.null(hit)) sym[names(hit)[!is.na(hit)]] <- hit[!is.na(hit)]
}
symbols <- coalesce(unname(sym[rownames(txi$counts)]), rownames(txi$counts))
names(symbols) <- rownames(txi$counts)
message("Genes without a symbol (shown as Ensembl ID): ",
        sum(symbols == names(symbols)), " / ", length(symbols))

# --- model -----------------------------------------------------------------
dge <- DGEList(round(txi$counts), samples = samples)
design <- model.matrix(as.formula(cfg$design$formula), data = samples)
dge <- dge[filterByExpr(dge, design), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
v <- voom(dge, design)
fit <- lmFit(v, design)
message("Model: ", nrow(dge), " genes x ", ncol(dge), " samples")

# Batch-corrected log-CPM — for visualisation only, never for testing.
logcpm <- if ("batch" %in% colnames(samples))
  removeBatchEffect(v$E, batch = samples$batch,
                    design = model.matrix(as.formula(paste("~", cond)), samples)) else v$E

saveRDS(list(fit = fit, design = design, dge = dge, logcpm = logcpm,
             samples = samples, symbols = symbols), obj_path("model.rds"))

# --- QC figures ------------------------------------------------------------
pca_plot <- function(mat, title) {
  sel <- order(rowVars(mat), decreasing = TRUE)[1:min(500, nrow(mat))]
  p <- prcomp(t(mat[sel, ]))
  pv <- round(100 * p$sdev^2 / sum(p$sdev^2))
  plot_data <- data.frame(
    PC1 = p$x[, 1], PC2 = p$x[, 2], samples[colnames(mat), ]) |>
    mutate(Condition = display_condition_factor(.data[[cond]]))
  ggplot(plot_data,
         aes(PC1, PC2, colour = Condition, label = colnames(mat))) +
    geom_point(size = 3) + ggrepel::geom_text_repel(size = 3, show.legend = FALSE) +
    scale_colour_manual(values = condition_palette(), drop = FALSE) + theme_bw() +
    labs(title = title, x = paste0("PC1 (", pv[1], "%)"),
         y = paste0("PC2 (", pv[2], "%)"), colour = "Condition")
}
save_plot(pca_plot(v$E, "PCA — before batch correction"), "pca_raw.svg")
save_plot(pca_plot(logcpm, "PCA — batch corrected"), "pca.svg")

d <- dist(t(logcpm))
distance_annotation <- data.frame(
  Condition = display_condition_factor(samples[[cond]]),
  row.names = rownames(samples))
save_plot(pheatmap::pheatmap(as.matrix(d), clustering_distance_rows = d,
  clustering_distance_cols = d, annotation_col = distance_annotation,
  annotation_colors = list(Condition = condition_palette()),
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(100),
  main = "Sample distances", silent = TRUE), "sample_distances.svg", 7, 6)

# Record the R platform and library paths used to create this run.
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))

message("01 done -> ", OUT)
