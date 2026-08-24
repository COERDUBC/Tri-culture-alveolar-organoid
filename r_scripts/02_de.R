#!/usr/bin/env Rscript
# 02 — Differential expression: one contrast per comparison in config.yml.
# Batch is already a covariate in the fitted model, so contrasts are batch-adjusted.
# Inputs: Stage 01's limma-voom model object plus comparison/threshold settings
# from config.yml. Outputs: all-gene DE tables, a significant-gene summary table,
# per-contrast volcano/heatmap figures, and targeted marker expression/effect plots.

source("r_scripts/common.R")
suppressPackageStartupMessages({ library(limma); library(pheatmap); library(tidyr) })

# model.rds contains the voom expression matrix, design matrix, fitted limma model,
# sample metadata, and Ensembl-to-symbol lookup created by 01_import.R.
m <- readRDS(obj_path("model.rds"))
cond <- cfg$design$condition_column
pad <- cfg$thresholds$padj; lfc <- cfg$thresholds$lfc

de_list <- lapply(cfg$comparisons, function(cmp) {
  # Build the configured treatment contrast. If the denominator is the model
  # reference level, limma already has a single coefficient for numerator vs
  # reference; otherwise makeContrasts subtracts the two condition coefficients.
  expr <- if (cmp$denominator == cfg$design$reference_level)
    paste0(cond, cmp$numerator) else
    paste0(cond, cmp$numerator, " - ", cond, cmp$denominator)

  # Apply the contrast to the stage-01 linear model and use empirical Bayes
  # moderation to stabilize gene-wise variance estimates across all genes.
  cv <- makeContrasts(contrasts = expr, levels = m$design)
  contrast_fit <- eBayes(contrasts.fit(m$fit, cv))

  # Convert limma's moderated variance pieces into a standard error and 95% CI
  # for the log2 fold-change estimate. These intervals are used for marker plots.
  moderated_se <- drop(contrast_fit$stdev.unscaled) *
    sqrt(contrast_fit$s2.post)
  uncertainty <- tibble(
    gene_id = rownames(contrast_fit$coefficients),
    std_error = moderated_se,
    residual_df = contrast_fit$df.total,
    conf_low = drop(contrast_fit$coefficients) -
      qt(0.975, contrast_fit$df.total) * moderated_se,
    conf_high = drop(contrast_fit$coefficients) +
      qt(0.975, contrast_fit$df.total) * moderated_se)

  # topTable reports one row per tested gene. log2FC is positive when expression
  # is higher in the numerator condition; adjusted P values use limma's default
  # multiple-testing correction because adjust.method is not set here.
  tt <- topTable(contrast_fit, number = Inf, sort.by = "P") |>
    rownames_to_column("gene_id") |> as_tibble() |>
    transmute(gene_id,
              gene = coalesce(unname(m$symbols[gene_id]), gene_id),
              log2FC = logFC, mean_logCPM = AveExpr, t = t,
              pvalue = P.Value, padj = adj.P.Val,
              significant = !is.na(padj) & padj < pad & abs(log2FC) > lfc,
              direction = case_when(significant & log2FC > 0 ~ "Up",
                                    significant & log2FC < 0 ~ "Down",
                                    TRUE ~ "Not significant")) |>
    left_join(uncertainty, by = "gene_id")
  write_csv(tt, tbl_path(paste0("de_", cmp$name, ".csv")))
  message(cmp$name, ": ", sum(tt$significant), " significant (padj<", pad,
          ", |log2FC|>", lfc, ")")
  tt
})
names(de_list) <- vapply(cfg$comparisons, function(x) x$name, character(1))
saveRDS(de_list, obj_path("de_list.rds"))

# One supplementary-ready table of every significant gene across the three
# contrasts. The per-contrast de_*.csv files above keep all tested genes; this is
# the same rows filtered to significant == TRUE, so no statistic is recomputed.
significant_de <- bind_rows(lapply(names(de_list), function(name)
  de_list[[name]] |>
    filter(significant) |>
    mutate(comparison = name,
           comparison_label = display_comparison(name), .before = 1))) |>
  arrange(comparison, padj, desc(abs(log2FC))) |>
  select(comparison, comparison_label, gene, gene_id, log2FC, conf_low,
         conf_high, mean_logCPM, t, pvalue, padj, direction)
write_csv(significant_de, tbl_path("de_significant_all_contrasts.csv"))
message("Significant genes across contrasts: ", nrow(significant_de))

# --- volcano + heatmap per comparison --------------------------------------
for (i in seq_along(de_list)) {
  name <- names(de_list)[i]
  d <- de_list[[name]]

  # Genes flagged in config are always labelled, significant or not (e.g. CXCL8
  # for Bleomycin, to place the transcript next to the measured IL-8 protein level).
  want <- cfg$comparisons[[i]]$label_genes %||% character()
  hl <- d |> filter(gene %in% want)
  if (length(setdiff(want, hl$gene)))
    warning(name, ": label_genes not in the model: ",
            paste(setdiff(want, hl$gene), collapse = ", "), call. = FALSE)
  if (nrow(hl))
    message(name, " labelled: ",
            paste0(hl$gene, " (log2FC=", round(hl$log2FC, 2),
                   ", padj=", signif(hl$padj, 2), ")", collapse = "; "))
  lab <- d |> filter(significant) |> slice_min(padj, n = 15) |>
    anti_join(hl, by = "gene_id")

  # Volcano plot: x-axis is effect size, y-axis is multiple-testing-adjusted
  # evidence. Dashed lines show the configured significance and fold-change cutoffs.
  p <- ggplot(d, aes(log2FC, -log10(pmax(padj, 1e-300)), colour = direction)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_hline(yintercept = -log10(pad), linetype = "dashed", linewidth = 0.3) +
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", linewidth = 0.3) +
    geom_point(data = hl, shape = 21, size = 2.8, stroke = 0.7,
               colour = "black", fill = NA, show.legend = FALSE) +
    # min.segment.length = 0 forces a leader line to every label, and the padding
    # plus force values push labels off the dense point cloud so the gene name and
    # the point it belongs to stay readable. seed keeps the layout reproducible.
    ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 3,
                             max.overlaps = Inf, show.legend = FALSE,
                             min.segment.length = 0, segment.size = 0.3,
                             segment.colour = "grey40", segment.alpha = 0.8,
                             box.padding = 0.7, point.padding = 0.35,
                             force = 4, force_pull = 0.4, seed = 42) +
    ggrepel::geom_text_repel(data = hl, aes(label = gene), size = 3.4,
                             fontface = "bold", colour = "black",
                             max.overlaps = Inf, show.legend = FALSE,
                             min.segment.length = 0, segment.size = 0.4,
                             segment.colour = "black",
                             box.padding = 1.0, point.padding = 0.5,
                             force = 6, seed = 42) +
    scale_colour_manual(values = c(Up = "#E64B35", Down = "#3C5488",
                                   `Not significant` = "grey70")) +
    labs(title = display_comparison(name), x = "log2 fold change",
         y = "-log10 adjusted p", colour = NULL) +
    theme_bw() + theme(legend.position = "bottom")
  save_plot(p, paste0("volcano_", name, ".svg"))

  # Heatmap: show the top significant genes for this contrast using the
  # batch-adjusted logCPM matrix saved for visualization, then z-score each gene
  # across samples so rows show relative within-gene patterns.
  top <- d |> filter(significant) |> slice_min(padj, n = cfg$thresholds$top_n_heatmap) |>
    pull(gene_id) |> intersect(rownames(m$logcpm))
  if (length(top) >= 2) {
    mat <- t(scale(t(m$logcpm[top, , drop = FALSE])))
    mat <- mat[rowSums(is.finite(mat)) == ncol(mat), , drop = FALSE]
    heatmap_annotation <- data.frame(
      Condition = display_condition_factor(m$samples[[cond]]),
      row.names = rownames(m$samples))
    save_plot(pheatmap(mat, annotation_col = heatmap_annotation,
      annotation_colors = list(Condition = condition_palette()),
      labels_row = coalesce(unname(m$symbols[rownames(mat)]), rownames(mat)),
      color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
      fontsize_row = 7, main = paste("Top genes:", display_comparison(name)),
      silent = TRUE), paste0("heatmap_", name, ".svg"), 7, 9)
  }
}

# --- Targeted aberrant-basal/alveolar marker figures -----------------------
target_genes <- as.character(unlist(cfg$targeted_marker_figure$genes))

# Map requested marker symbols back to model rows. When multiple Ensembl rows map
# to the same symbol, keep the row with the highest mean plotted expression.
marker_map <- tibble(
  gene_id = rownames(m$logcpm),
  gene = coalesce(unname(m$symbols[rownames(m$logcpm)]), rownames(m$logcpm)),
  mean_log2_expression = rowMeans(m$logcpm)) |>
  filter(gene %in% target_genes) |>
  arrange(match(gene, target_genes), desc(mean_log2_expression), gene_id) |>
  distinct(gene, .keep_all = TRUE)
missing_markers <- setdiff(target_genes, marker_map$gene)
# A missing gene gets no facet, so the warning is the only record of it. SFTPC is
# the expected case — the AT2 component is the NCI-H441 line, which does not
# express mature surfactant protein C.
if (length(missing_markers))
  warning("Targeted marker genes absent after expression filtering: ",
          paste(missing_markers, collapse = ", "), call. = FALSE)
if (nrow(marker_map) < 2L)
  stop("Fewer than two targeted marker genes remain in the model.",
       call. = FALSE)
write_csv(marker_map, tbl_path("targeted_marker_gene_mapping.csv"))

# Long-format per-sample marker expression for plotting individual sample points
# and condition summaries. These are descriptive values, not a second DE test.
marker_expression <- as.data.frame(
  m$logcpm[marker_map$gene_id, , drop = FALSE], check.names = FALSE) |>
  rownames_to_column("gene_id") |>
  pivot_longer(-gene_id, names_to = "sample", values_to = "log2_expression") |>
  left_join(select(marker_map, gene_id, gene), by = "gene_id") |>
  mutate(
    condition = as.character(m$samples[sample, cond]),
    Condition = display_condition_factor(condition),
    gene = factor(gene, levels = target_genes))
marker_summary <- marker_expression |>
  group_by(gene, condition, Condition) |>
  summarise(
    mean_log2_expression = mean(log2_expression),
    sem = sd(log2_expression) / sqrt(n()), n = n(), .groups = "drop")
write_csv(marker_expression, tbl_path("targeted_marker_expression.csv"))
write_csv(marker_summary,
          tbl_path("targeted_marker_expression_summary.csv"))

# Reuse the same contrast-level DE results for targeted genes, so the marker
# effect plot and significance labels stay aligned with the all-gene analysis.
marker_effects <- bind_rows(lapply(names(de_list), function(name)
  de_list[[name]] |>
    semi_join(select(marker_map, gene_id), by = "gene_id") |>
    mutate(comparison = name, .before = 1))) |>
  mutate(
    comparison_label = factor(
      display_comparison(comparison),
      levels = display_comparison(names(de_list))),
    gene = factor(gene, levels = target_genes),
    significance = case_when(
      padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*",
      TRUE ~ "ns"))
write_csv(marker_effects, tbl_path("targeted_marker_treatment_effects.csv"))

# Add all three configured pairwise annotations to the descriptive expression
# figure. P values come from the same limma-voom result object written above;
# no second test is performed on the batch-adjusted expression values.
comparison_key <- function(x)
  paste(as.character(x$numerator), as.character(x$denominator), sep = "_vs_")
# Nest the two adjacent comparisons below the longest Control-to-Bleomycin
# bracket. These are display layers only; all P values come from the saved
# limma-voom treatment-effects table above.
annotation_order <- c("TGF_vs_Control", "TGF_vs_Bleo", "Bleo_vs_Control")
configured_keys <- vapply(cfg$comparisons, comparison_key, character(1))
expression_comparison_cfg <- cfg$comparisons[
  match(annotation_order, configured_keys)]
if (length(expression_comparison_cfg) != 3L ||
    any(vapply(expression_comparison_cfg, is.null, logical(1))))
  stop("Targeted expression annotations require TGF-vs-Control, ",
       "TGF-vs-Bleomycin, and Bleomycin-vs-Control comparisons in config.yml.",
       call. = FALSE)
condition_levels <- levels(marker_expression$Condition)

# Translate configured denominator/numerator pairs into x-axis positions on the
# expression plot. The bracket layer controls vertical stacking within each facet.
annotation_spec <- bind_rows(lapply(seq_along(expression_comparison_cfg),
  function(i) {
    cmp <- expression_comparison_cfg[[i]]
    tibble(
      comparison = as.character(cmp$name),
      numerator = as.character(cmp$numerator),
      denominator = as.character(cmp$denominator),
      bracket_layer = i,
      x_start = min(match(display_condition(cmp$denominator), condition_levels),
                    match(display_condition(cmp$numerator), condition_levels)),
      x_end = max(match(display_condition(cmp$denominator), condition_levels),
                  match(display_condition(cmp$numerator), condition_levels)))
  }))
if (anyNA(annotation_spec$x_start) || anyNA(annotation_spec$x_end))
  stop("Could not map targeted-expression comparison labels to plot positions.",
       call. = FALSE)

expected_annotations <- tidyr::expand_grid(
  gene = as.character(marker_map$gene),
  comparison = annotation_spec$comparison)
missing_annotations <- anti_join(
  expected_annotations,
  filter(marker_effects, is.finite(padj)) |>
    distinct(gene, comparison),
  by = c("gene", "comparison"))
if (nrow(missing_annotations))
  warning("Targeted expression significance omitted for missing/non-finite ",
          "adjusted P value(s): ",
          paste(paste(missing_annotations$gene,
                      missing_annotations$comparison, sep = "/"),
                collapse = ", "), call. = FALSE)

expression_extents <- marker_expression |>
  mutate(gene = as.character(gene)) |>
  group_by(gene) |>
  summarise(
    y_min = min(log2_expression, na.rm = TRUE),
    sample_max = max(log2_expression, na.rm = TRUE),
    .groups = "drop") |>
  left_join(
    marker_summary |>
      mutate(gene = as.character(gene)) |>
      group_by(gene) |>
      summarise(
        summary_max = max(mean_log2_expression + coalesce(sem, 0),
                          na.rm = TRUE),
        .groups = "drop"),
    by = "gene") |>
  mutate(
    y_max = pmax(sample_max, summary_max),
    # Bars start at 0, so a panel spans [min(0, y_min), y_max] — offsets have to be
    # a fraction of THAT, not of the data spread. Scaling by the spread gives panels
    # whose values sit far from zero (CDKN1A, SLC2A1: spread ~1 in a panel ~10 tall)
    # brackets a few pixels apart that clip out of the panel, while panels reaching
    # down toward zero (KRT17) get roomy ones. pmin(0, .) is for MMP7, which has a
    # sample below zero.
    panel_bottom = pmin(0, y_min),
    panel_span = pmax(y_max - panel_bottom, 0.5),
    bracket_base = y_max + 0.06 * panel_span,
    bracket_step = 0.09 * panel_span,
    bracket_tick = 0.02 * panel_span)

# Attach each marker's adjusted P value to its plotting coordinates and convert
# FDR thresholds to star labels. Non-significant comparisons are filtered below.
expression_annotations <- marker_effects |>
  filter(comparison %in% annotation_spec$comparison, is.finite(padj)) |>
  select(gene, comparison, padj) |>
  inner_join(annotation_spec, by = "comparison") |>
  inner_join(expression_extents, by = "gene") |>
  mutate(
    gene = factor(gene, levels = target_genes),
    significance = case_when(
      padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*",
      TRUE ~ "ns"),
    bracket_y = bracket_base + (bracket_layer - 1) * bracket_step,
    tick_y = bracket_y - bracket_tick,
    label_x = (x_start + x_end) / 2)

# Only significant comparisons are drawn. A bare "ns" with no bracket has
# nothing to attach it to, so the whole annotation goes.
expression_brackets <- filter(expression_annotations, significance != "ns")

# Descriptive marker expression figure: bars are condition means, error bars are
# SEM, points are samples, and significance brackets are imported from limma DE.
save_plot(
  ggplot(marker_summary,
         aes(Condition, mean_log2_expression, fill = Condition)) +
    geom_col(width = 0.68, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = mean_log2_expression - sem,
          ymax = mean_log2_expression + sem), width = 0.18) +
    geom_point(
      data = marker_expression,
      aes(Condition, log2_expression), inherit.aes = FALSE,
      position = position_jitter(width = 0.08, height = 0, seed = 20260715),
      shape = 21, size = 2, stroke = 0.45, fill = "white") +
    geom_segment(
      data = expression_brackets,
      aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(
      data = expression_brackets,
      aes(x = x_start, xend = x_start, y = bracket_y, yend = tick_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(
      data = expression_brackets,
      aes(x = x_end, xend = x_end, y = bracket_y, yend = tick_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    # vjust, not a data-unit offset: each facet has its own free y scale, so a
    # fixed offset in expression units lands at a different distance in every
    # panel. vjust is in text heights, which is identical everywhere.
    geom_text(
      data = expression_brackets,
      aes(x = label_x, y = bracket_y, label = significance),
      inherit.aes = FALSE, size = 3, vjust = -0.15) +
    facet_wrap(~gene, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = condition_palette(), drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    labs(
      title = "Fibrosis and alveolar marker expression",
      x = NULL, y = "Batch-adjusted voom log2 expression",
      fill = "Condition") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(face = "italic"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom"),
  "targeted_aberrant_basal_marker_expression.svg", 12, 7)

comparison_colours <- setNames(
  PAL[seq_along(levels(marker_effects$comparison_label))],
  levels(marker_effects$comparison_label))

# Targeted marker effect figure: show the model-estimated log2FC and 95% CI for
# each configured contrast, rather than raw group means.
save_plot(
  ggplot(marker_effects,
         aes(log2FC, comparison_label, colour = comparison_label)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey55") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high),
                  width = 0.18, orientation = "y") +
    geom_point(size = 2.4) +
    geom_text(aes(x = conf_high, label = significance),
              hjust = -0.35, colour = "black", size = 3) +
    facet_wrap(~gene, ncol = 5) +
    scale_colour_manual(values = comparison_colours, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.20))) +
    labs(
      title = "Targeted marker contrasts",
      x = "log2 fold change (numerator − denominator)", y = NULL,
      colour = "Comparison") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(face = "italic"),
          legend.position = "bottom"),
  "targeted_aberrant_basal_marker_effects.svg", 12, 6.5)

message("02 done -> ", OUT)
