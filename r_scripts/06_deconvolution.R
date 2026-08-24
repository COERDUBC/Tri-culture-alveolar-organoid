#!/usr/bin/env Rscript
# 06 — Lineage/reference contributions in the tri-culture.
#
# Lineage-restricted GSE136831 CIBERSORTx (five states), aggregated to
# AT1 / Epithelial / Fibroblast for the main comparison and downstream use.
# A second reference built with a smaller donor cap is read as a sensitivity.
#
# These are RNA/reference contribution scores, not validated cell counts.

source("r_scripts/common.R")
suppressPackageStartupMessages({
  library(edgeR); library(limma); library(tidyr)
})

m <- readRDS(obj_path("model.rds"))
cond <- cfg$design$condition_column

bulk <- cpm(to_symbol(m$dge$counts), log = FALSE)
expected_samples <- colnames(bulk)

required_file <- function(path, label) {
  missing <- is.null(path) || length(path) != 1 || is.na(path) ||
    !nzchar(trimws(path)) || grepl("^TO_FILL$", trimws(path)) ||
    !file.exists(path) || dir.exists(path)
  if (missing) stop(label, " is unavailable: ", path %||% "<unset>", call. = FALSE)
  path
}

long <- function(x, method) as.data.frame(x) |>
  rownames_to_column("sample") |>
  pivot_longer(-sample, names_to = "cell_type", values_to = "contribution") |>
  mutate(method = method, .before = 1)

read_cibersortx <- function(path, label) {
  raw <- read_csv(path, show_col_types = FALSE)
  sample_col <- intersect(c("Mixture", "Sample", "sample"), names(raw))
  if (!length(sample_col)) sample_col <- names(raw)[1]
  samples <- as.character(raw[[sample_col[1]]])
  if (anyNA(samples) || any(samples == "") || anyDuplicated(samples))
    stop(label, " has blank or duplicated sample identifiers.", call. = FALSE)
  if (!setequal(samples, expected_samples))
    stop(label, " samples do not match bulk samples. ",
         set_diff_text(expected_samples, samples), call. = FALSE)

  # CIBERSORTx labels the absolute-mode total differently across exports
  # (e.g. "Absolute score" or "Absolute score (sig.score)").
  diagnostic_cols <- names(raw)[grepl(
    "^(P[- .]?value|Correlation|RMSE|Absolute[ .]score)([[:space:]]*\\([^)]*\\))?$",
    names(raw), ignore.case = TRUE, perl = TRUE)]
  
  fraction_cols <- setdiff(names(raw), c(sample_col[1], diagnostic_cols))
  fractions <- as.matrix(data.frame(
    lapply(raw[fraction_cols], function(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE))
  rownames(fractions) <- samples
  fractions <- fractions[expected_samples, , drop = FALSE]
  if (!ncol(fractions) || any(!is.finite(fractions)) ||
      any(fractions < -1e-8))
    stop(label, " has non-numeric, non-finite, or negative fractions.",
         call. = FALSE)

  diagnostics <- tibble(sample = samples)
  if (length(diagnostic_cols))
    diagnostics <- bind_cols(
      diagnostics,
      raw[, diagnostic_cols, drop = FALSE])
  diagnostics <- diagnostics[match(expected_samples, diagnostics$sample), ,
                             drop = FALSE]
  list(fractions = fractions, diagnostics = diagnostics, raw = raw)
}

resolve_groups <- function(cell_types, patterns) {
  groups <- lapply(patterns, function(pattern)
    cell_types[grepl(pattern, cell_types, ignore.case = TRUE, perl = TRUE)])
  if (any(lengths(groups) == 0))
    stop("No GSE136831 state matched compartment(s): ",
         paste(names(groups)[lengths(groups) == 0], collapse = ", "),
         call. = FALSE)
  assigned <- unlist(groups, use.names = FALSE)
  overlap <- unique(assigned[duplicated(assigned)])
  unassigned <- setdiff(cell_types, assigned)
  if (length(overlap))
    stop("GSE136831 states matched multiple compartments: ",
         paste(overlap, collapse = ", "), call. = FALSE)
  if (length(unassigned))
    stop("GSE136831 states were not assigned to a broad compartment: ",
         paste(unassigned, collapse = ", "), call. = FALSE)
  groups
}

# --- Required CIBERSORTx results -------------------------------------------
cib_cfg <- cfg$deconvolution$cibersortx
GSE136831_path <- required_file(
  cib_cfg$GSE136831_lineage_result,
  "GSE136831 lineage CIBERSORTx result")
GSE136831_cib <- read_cibersortx(
  GSE136831_path, "GSE136831 lineage CIBERSORTx result")

# CIBERSORTx absolute mode scales each sample by its estimated total mRNA
# content, so rows do not sum to 1. Renormalizing would convert the scores back
# into compositions and discard exactly what absolute mode provides, so both
# scales are kept: absolute for comparing a cell type ACROSS samples, relative
# for describing the composition WITHIN a sample.
absolute <- isTRUE(cfg$deconvolution$cibersortx$absolute_mode)
message("CIBERSORTx mode: ", if (absolute) "absolute" else "relative")

expected_states <- unique(as.character(unlist(
  cfg$deconvolution$GSE136831_atlas$lineage_restricted_cell_types)))
if (!setequal(colnames(GSE136831_cib$fractions), expected_states))
  stop("GSE136831 CIBERSORTx columns do not match configured lineage states. ",
       set_diff_text(expected_states, colnames(GSE136831_cib$fractions)),
       call. = FALSE)
GSE136831_states <- GSE136831_cib$fractions[, expected_states, drop = FALSE]
GSE136831_state_sums <- rowSums(GSE136831_states)
if (any(GSE136831_state_sums <= 0))
  stop("GSE136831 five-state scores contain a zero-total sample.",
       call. = FALSE)
GSE136831_state_rel <- GSE136831_states / GSE136831_state_sums

groups <- resolve_groups(
  expected_states,
  cfg$deconvolution$GSE136831_atlas$comparison_patterns)
GSE136831_three_raw <- sapply(groups, function(states)
  rowSums(GSE136831_states[, states, drop = FALSE]))
rownames(GSE136831_three_raw) <- expected_samples
selected_mass <- rowSums(GSE136831_three_raw)
total_mass <- rowSums(GSE136831_states)
if (any(selected_mass <= 0) || any(!is.finite(selected_mass)))
  stop("A sample has zero/non-finite GSE136831 broad-compartment mass.",
       call. = FALSE)
compartments <- c("AT1", "Epithelial", "Fibroblast")
GSE136831_abs <- GSE136831_three_raw[, compartments, drop = FALSE]
GSE136831_rel <- GSE136831_abs / selected_mass
GSE136831_three <- if (absolute) GSE136831_abs else GSE136831_rel
if (max(abs(selected_mass - total_mass)) > 1e-8)
  stop("Broad-compartment aggregation did not retain every GSE136831 state.",
       call. = FALSE)
# Row sums are only expected to be 1 in relative mode. In absolute mode the row
# sum IS the CIBERSORTx absolute score, so it is recorded rather than flagged.
if (!absolute && max(abs(total_mass - 1)) > 0.02)
  warning("GSE136831 CIBERSORTx row sums differ from 1 by more than 0.02 in ",
          "relative mode; raw sums are retained in QC. If the job was run in ",
          "absolute mode, set deconvolution$cibersortx$absolute_mode: true.",
          call. = FALSE)
if (absolute)
  message("Absolute scores (row sums): ",
          paste(round(range(total_mass), 3), collapse = " - "))

# --- Batch-adjusted treatment comparisons ----------------------------------
# Absolute-mode scores are not constrained to sum to one, so each state or
# broad compartment can be compared across conditions without inducing a
# compositional decrease in every other component. Each score receives its own
# ordinary linear model; no empirical-Bayes borrowing is performed across the
# small, biologically dependent set of reference components.
test_score_matrix <- function(scores, resolution) {
  if (!identical(rownames(scores), expected_samples))
    stop("Treatment-test score rows are not in bulk-sample order.",
         call. = FALSE)
  sample_data <- m$samples[expected_samples, , drop = FALSE]
  design <- model.matrix(as.formula(cfg$design$formula), data = sample_data)
  if (qr(design)$rank < ncol(design))
    stop("Deconvolution treatment model is rank deficient.", call. = FALSE)
  fit <- lmFit(t(scores), design)

  bind_rows(lapply(cfg$comparisons, function(cmp) {
    numerator <- as.character(cmp$numerator)
    denominator <- as.character(cmp$denominator)
    expression <- if (denominator == cfg$design$reference_level)
      paste0(cond, numerator) else
      paste0(cond, numerator, " - ", cond, denominator)
    tested <- contrasts.fit(
      fit, makeContrasts(contrasts = expression, levels = design))
    estimate <- drop(tested$coefficients)
    std_error <- drop(tested$stdev.unscaled) * tested$sigma
    residual_df <- tested$df.residual
    statistic <- estimate / std_error
    margin <- qt(0.975, residual_df) * std_error

    tibble(
      resolution, cell_type = colnames(scores), contrast = cmp$name,
      mean_numerator = colMeans(
        scores[sample_data[[cond]] == numerator, , drop = FALSE]),
      mean_denominator = colMeans(
        scores[sample_data[[cond]] == denominator, , drop = FALSE]),
      adjusted_difference = estimate, std_error,
      conf_low = estimate - margin, conf_high = estimate + margin,
      statistic, residual_df,
      pvalue = 2 * pt(-abs(statistic), residual_df)
    )
  })) |>
    group_by(resolution, contrast) |>
    mutate(
      # 0.05 is deliberately not cfg$thresholds$padj: that key is the DE gate
      # (Stage 02/03). This is the deconvolution score-model gate, and the two
      # are kept free to move independently.
      padj = p.adjust(pvalue, method = "BH"),
      significant = !is.na(padj) & padj < 0.05
    ) |>
    ungroup()
}

if (!absolute)
  stop("Batch-adjusted cross-treatment deconvolution tests require ",
       "deconvolution$cibersortx$absolute_mode: true. Relative fractions are ",
       "compositional and are retained for within-sample description only.",
       call. = FALSE)

treatment_tests <- bind_rows(
  test_score_matrix(GSE136831_states, "five_state"),
  test_score_matrix(GSE136831_abs, "three_compartment")
) |>
  arrange(resolution, contrast, padj, cell_type)

# --- Reference-size sensitivity --------------------------------------------
# A second CIBERSORTx result built from a larger single-cell reference (more
# cells per donor). It is pushed through the SAME test_score_matrix() call, so
# any difference is attributable to the reference and not to the analysis.
# Answers the review question of whether the treatment effects survive a bigger
# reference; the primary result above is unaffected either way.
sensitivity_path <- cib_cfg$GSE136831_lineage_result_sensitivity %||% ""
reference_sensitivity <- NULL
if (nzchar(trimws(sensitivity_path)) &&
    !grepl("^TO_FILL$", trimws(sensitivity_path))) {
  sens_cib <- read_cibersortx(
    required_file(sensitivity_path, "GSE136831 sensitivity CIBERSORTx result"),
    "GSE136831 sensitivity CIBERSORTx result")
  if (!setequal(colnames(sens_cib$fractions), expected_states))
    stop("Sensitivity CIBERSORTx columns do not match configured lineage ",
         "states. ",
         set_diff_text(expected_states, colnames(sens_cib$fractions)),
         call. = FALSE)
  sens_states <- sens_cib$fractions[, expected_states, drop = FALSE]

  reference_sensitivity <- treatment_tests |>
    filter(resolution == "five_state") |>
    select(cell_type, contrast,
           primary_difference = adjusted_difference,
           primary_padj = padj, primary_significant = significant) |>
    inner_join(
      test_score_matrix(sens_states, "five_state") |>
        select(cell_type, contrast,
               sensitivity_difference = adjusted_difference,
               sensitivity_padj = padj,
               sensitivity_significant = significant),
      by = c("cell_type", "contrast")) |>
    mutate(
      same_direction = sign(primary_difference) == sign(sensitivity_difference),
      same_significance = primary_significant == sensitivity_significant,
      agree = same_direction & same_significance) |>
    arrange(contrast, cell_type)
  write_csv(reference_sensitivity,
            tbl_path("deconvolution_reference_sensitivity.csv"))

  # Separates "the absolute scores moved" from "the conclusions moved": a
  # reference change can rescale every score without altering any contrast.
  write_csv(
    tibble(
      cell_type = expected_states,
      pearson = vapply(expected_states, function(state)
        cor(GSE136831_states[, state], sens_states[, state]), numeric(1)),
      spearman = vapply(expected_states, function(state)
        cor(GSE136831_states[, state], sens_states[, state],
            method = "spearman"), numeric(1)),
      primary_mean = colMeans(GSE136831_states)[expected_states],
      sensitivity_mean = colMeans(sens_states)[expected_states]),
    tbl_path("deconvolution_reference_score_agreement.csv"))

  # CIBERSORTx's own per-sample goodness of fit for each reference. This is the
  # criterion for choosing a primary, because it never looks at condition.
  write_csv(
    bind_rows(
      mutate(GSE136831_cib$diagnostics, reference = "primary", .before = 1),
      mutate(sens_cib$diagnostics, reference = "sensitivity", .before = 1)),
    tbl_path("deconvolution_reference_fit_quality.csv"))

  message("Reference sensitivity: ", sum(reference_sensitivity$agree), " / ",
          nrow(reference_sensitivity),
          " state x contrast results agree in direction and significance.")
}

# --- Outputs ---------------------------------------------------------------
# `scale` distinguishes the two views: absolute is comparable across samples,
# relative describes the composition within a sample. The primary scale is the
# one matching the CIBERSORTx mode; both are written so either can be plotted.
primary_scale <- if (absolute) "absolute" else "relative"
primary_long <- long(GSE136831_states, "CIBERSORTx_GSE136831_5state_primary") |>
  mutate(scale = primary_scale, .after = method) |>
  mutate(!!cond := m$samples[sample, cond],
         Condition = display_condition_factor(.data[[cond]]))
both_scales_long <- bind_rows(
  long(GSE136831_states, "CIBERSORTx_GSE136831_5state_primary") |>
    mutate(scale = "absolute"),
  long(GSE136831_state_rel, "CIBERSORTx_GSE136831_5state_primary") |>
    mutate(scale = "relative")) |>
  relocate(scale, .after = method) |>
  mutate(!!cond := m$samples[sample, cond])
state_long <- long(GSE136831_states, "CIBERSORTx_GSE136831_state") |>
  mutate(!!cond := m$samples[sample, cond],
         Condition = display_condition_factor(.data[[cond]]),
         cell_type = factor(cell_type, levels = expected_states))
write_csv(both_scales_long, tbl_path("deconvolution_both_scales.csv"))
write_csv(state_long, tbl_path("deconvolution_GSE136831_states.csv"))
write_csv(
  long(GSE136831_three, "CIBERSORTx_GSE136831_3compartment_QC") |>
    mutate(!!cond := m$samples[sample, cond]),
  tbl_path("deconvolution_three_compartment_qc.csv"))
write_csv(as.data.frame(GSE136831_three_raw) |>
            rownames_to_column("sample"),
          tbl_path("deconvolution_GSE136831_three_raw.csv"))
write_csv(tibble(
  sample = expected_samples,
  cibersortx_mode = primary_scale,
  # In absolute mode this row sum IS the CIBERSORTx absolute score.
  selected_three_group_mass = selected_mass,
  total_lineage_mass = total_mass,
  aggregation_difference = selected_mass - total_mass),
  tbl_path("deconvolution_GSE136831_coverage.csv"))
write_csv(tibble(
  compartment = names(groups),
  lineage_states = vapply(groups, paste, character(1), collapse = ";")),
  tbl_path("deconvolution_GSE136831_group_mapping.csv"))
write_csv(GSE136831_cib$diagnostics,
          tbl_path("deconvolution_GSE136831_diagnostics.csv"))
write_csv(treatment_tests,
          tbl_path("deconvolution_treatment_tests.csv"))

# Figure 1-style grouped summary. The plot is analogous to Ryu et al. rather
# than a reproduction: it shows the five prespecified organoid-relevant states
# and the three experimental conditions. Two bracket layers reuse the
# batch-adjusted score-model P values above; no test is run on plotted values.
deconv_comparison_cfg <- Filter(function(x)
  identical(as.character(x$denominator), cfg$design$reference_level) &&
    as.character(x$numerator) %in% c("TGF", "Bleo"), cfg$comparisons)
deconv_comparison_cfg <- deconv_comparison_cfg[
  order(match(vapply(deconv_comparison_cfg, function(x) x$numerator,
                     character(1)), c("TGF", "Bleo")))]
if (length(deconv_comparison_cfg) != 2L)
  stop("Five-state brackets require TGF-vs-Control and ",
       "Bleomycin-vs-Control comparisons in config.yml.", call. = FALSE)

dodge_width <- 0.76
condition_levels <- levels(state_long$Condition)
condition_offsets <- setNames(
  dodge_width * ((seq_along(condition_levels) - 0.5) /
                   length(condition_levels) - 0.5),
  condition_levels)
deconv_bracket_spec <- bind_rows(lapply(seq_along(deconv_comparison_cfg),
  function(i) {
    cmp <- deconv_comparison_cfg[[i]]
    tibble(
      contrast = as.character(cmp$name),
      bracket_layer = i,
      start_condition = display_condition(cmp$denominator),
      end_condition = display_condition(cmp$numerator))
  })) |>
  mutate(
    start_offset = unname(condition_offsets[start_condition]),
    end_offset = unname(condition_offsets[end_condition]))
if (anyNA(deconv_bracket_spec$start_offset) ||
    anyNA(deconv_bracket_spec$end_offset))
  stop("Could not align five-state brackets to condition boxplots.",
       call. = FALSE)

state_extents <- state_long |>
  mutate(cell_type = as.character(cell_type)) |>
  group_by(cell_type) |>
  summarise(data_max = max(contribution, na.rm = TRUE), .groups = "drop") |>
  mutate(
    state_x = match(cell_type, expected_states),
    bracket_base = pmin(data_max + 0.07, 1.28))

deconv_brackets <- treatment_tests |>
  filter(resolution == "five_state",
         contrast %in% deconv_bracket_spec$contrast) |>
  select(cell_type, contrast, padj) |>
  inner_join(deconv_bracket_spec, by = "contrast") |>
  inner_join(state_extents, by = "cell_type") |>
  mutate(
    cell_type = factor(cell_type, levels = expected_states),
    label = case_when(
      padj < 0.001 ~ "***", padj < 0.01 ~ "**", padj < 0.05 ~ "*",
      TRUE ~ "ns"),
    x_start = state_x + pmin(start_offset, end_offset),
    x_end = state_x + pmax(start_offset, end_offset),
    bracket_y = bracket_base + (bracket_layer - 1) * 0.12,
    tick_y = bracket_y - 0.025,
    label_x = (x_start + x_end) / 2,
    label_y = bracket_y + 0.025)
if (nrow(deconv_brackets) != length(expected_states) * 2L ||
    any(!is.finite(deconv_brackets$padj)))
  stop("Five-state bracket annotations require two finite adjusted P values ",
       "for every reference state.", call. = FALSE)
if (any(deconv_brackets$label_y > 1.5))
  stop("Five-state bracket annotations exceed the fixed 0-1.5 y-axis.",
       call. = FALSE)

save_plot(
  ggplot(state_long,
         aes(cell_type, contribution, fill = Condition)) +
    geom_vline(
      xintercept = seq(1.5, length(expected_states) - 0.5, by = 1),
      colour = "grey90", linewidth = 0.35) +
    geom_boxplot(
      position = position_dodge(width = dodge_width), width = 0.68,
      outlier.shape = NA, linewidth = 0.35) +
    geom_point(
      aes(colour = Condition),
      position = position_jitterdodge(
        jitter.width = 0.10, dodge.width = dodge_width),
      shape = 21, fill = "white", stroke = 0.35, size = 1.5,
      show.legend = FALSE) +
    geom_segment(
      data = deconv_brackets,
      aes(x = x_start, xend = x_end, y = bracket_y, yend = bracket_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(
      data = deconv_brackets,
      aes(x = x_start, xend = x_start, y = bracket_y, yend = tick_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    geom_segment(
      data = deconv_brackets,
      aes(x = x_end, xend = x_end, y = bracket_y, yend = tick_y),
      inherit.aes = FALSE, linewidth = 0.35) +
    geom_text(
      data = deconv_brackets,
      aes(x = label_x, y = label_y, label = label),
      inherit.aes = FALSE, size = 3) +
    scale_fill_manual(values = condition_palette(), drop = FALSE) +
    scale_colour_manual(values = condition_palette(), drop = FALSE) +
    scale_y_continuous(
      breaks = seq(0, 1.5, by = 0.5),
      expand = expansion(mult = c(0, 0.02))) +
    coord_cartesian(ylim = c(0, 1.5)) +
    theme_bw(base_size = 9) +
    labs(y = "CIBERSORTx absolute transcriptomic abundance score", x = NULL,
         title = "Five-state transcriptomic abundance scores",
         fill = "Condition") +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "right"),
  "deconvolution_GSE136831_states.svg", 9.5, 5.2)

effect_position <- position_dodge(width = 0.55)
save_plot(
  ggplot(filter(treatment_tests, resolution == "five_state") |>
           mutate(Comparison = factor(
             display_comparison(contrast),
             levels = display_comparison(vapply(
               cfg$comparisons, function(x) x$name, character(1))))),
         aes(adjusted_difference, cell_type, colour = Comparison,
             group = Comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3,
               colour = "grey45") +
    geom_errorbar(aes(xmin = conf_low, xmax = conf_high),
                  position = effect_position, width = 0.18,
                  orientation = "y") +
    geom_point(aes(shape = significant), position = effect_position, size = 2.6) +
    scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                       labels = c(`FALSE` = "BH-adjusted P ≥ 0.05",
                                  `TRUE` = "BH-adjusted P < 0.05")) +
    theme_bw(base_size = 9) +
    labs(
      title = "Batch-adjusted contrasts of five-state CIBERSORTx scores",
      x = "adjusted absolute-score difference (numerator − denominator)",
      y = NULL, colour = "Comparison", shape = NULL),
  "deconvolution_treatment_effects.svg", 9, 7)

# Archive one canonical field per result.
deconvolution_object <- list(
  primary_method = "CIBERSORTx_GSE136831_lineage",
  primary_five_state_scores = GSE136831_states,
  relative_reference_contributions = GSE136831_state_rel,
  cibersortx_GSE136831_lineage = GSE136831_cib,
  cibersortx_GSE136831_comparison_groups = GSE136831_three,
  cibersortx_GSE136831_comparison_groups_raw = GSE136831_three_raw,
  cibersortx_GSE136831_selected_mass = selected_mass,
  deconvolution_treatment_tests = treatment_tests,
  reference_size_sensitivity = reference_sensitivity,
  long = primary_long)
saveRDS(deconvolution_object, obj_path("deconvolution.rds"))

# boxplot for the deconvolution
plot_state_long <- state_long |>
  mutate(
    cell_type = factor(
      cell_type,
      levels = sort(expected_states)
    )
  )

save_plot(
  ggplot(
    plot_state_long,
    aes(x = Condition, y = contribution, fill = Condition)
  ) +
    geom_boxplot(
      outlier.size = 0.4,
      linewidth = 0.35
    ) +
    facet_wrap(
      ~cell_type,
      scales = "free_y",
      nrow = 1
    ) +
    scale_fill_manual(
      values = condition_palette(),
      drop = FALSE
    ) +
    theme_bw(base_size = 8) +
    labs(
      x = NULL,
      y = "CIBERSORTx absolute transcriptomic abundance score",
      fill = "Condition",
      title = paste(
        "Cellular deconvolution for five expected cell states",
        "in the organoids"
      )
    ) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ),
  "deconvolution_GSE136831_states_boxplot.svg",
  width = 12,
  height = 4.5
)

message("Primary: lineage-restricted GSE136831 CIBERSORTx.")
message("Batch-adjusted treatment tests: ",
        sum(treatment_tests$significant), " / ", nrow(treatment_tests),
        " rows have BH-adjusted P < 0.05.")
message("06 done -> ", OUT)
