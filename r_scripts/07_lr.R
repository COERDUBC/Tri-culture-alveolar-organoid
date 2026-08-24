#!/usr/bin/env Rscript
# 07 — Candidate cross-state ligand-receptor signalling (connectomeDB2025).
#
# Direction comes from the five-state GSE136831 single-cell reference. A gene
# may belong to multiple states when it passes both the mean-CPM and detection
# thresholds.

source("r_scripts/common.R")
suppressPackageStartupMessages({ library(edgeR); library(tidyr) })

m <- readRDS(obj_path("model.rds"))
cond <- cfg$design$condition_column
lr_csv <- cfg$crosstalk$lr_csv %||% ""
if (!nzchar(lr_csv) || !file.exists(lr_csv))
  stop("Set crosstalk$lr_csv to the connectomeDB2025 human CSV in config.yml.")

tri <- to_symbol(m$dge$counts)
cpm_lin <- cpm(tri, log = FALSE)
log_cpm <- cpm(tri, log = TRUE, prior.count = 0.5)
min_cpm <- cfg$crosstalk$lr_min_cpm %||% 1
expressed <- rownames(cpm_lin)[rowMeans(cpm_lin >= min_cpm) >= 0.5]

# --- input interactions ----------------------------------------------------
raw <- read_csv(lr_csv, show_col_types = FALSE)
lcol <- cfg$crosstalk$lr_ligand_col %||% "Ligand_Symbols"
rcol <- cfg$crosstalk$lr_receptor_col %||% "Receptor_Symbols"
if (!all(c(lcol, rcol) %in% names(raw)))
  stop("Columns ", lcol, "/", rcol, " not in ", basename(lr_csv), ": ",
       paste(names(raw), collapse = ", "))
approved <- function(x) sub("\\s*\\(.*$", "", trimws(as.character(x)))
lr_base <- tibble(ligand = approved(raw[[lcol]]),
                  receptor = approved(raw[[rcol]])) |>
  distinct() |>
  filter(ligand %in% expressed, receptor %in% expressed)

# Correlations are gene-pair properties and are calculated once before state
# expansion. Residual rho is the primary descriptive co-expression estimate.
resid <- t(lm.fit(model.matrix(as.formula(cfg$design$formula), m$samples),
                  t(log_cpm))$residuals)
dimnames(resid) <- dimnames(log_cpm)
condition <- as.character(m$samples[colnames(log_cpm), cond])
rho <- function(x, y) {
  if (length(x) < 4 || sd(x) == 0 || sd(y) == 0)
    return(c(NA_real_, NA_real_))
  z <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(unname(z$estimate), z$p.value)
}
pair_stats <- bind_rows(lapply(seq_len(nrow(lr_base)), function(i) {
  l <- log_cpm[lr_base$ligand[i], ]
  r <- log_cpm[lr_base$receptor[i], ]
  pooled <- rho(l, r)
  adjusted <- rho(resid[lr_base$ligand[i], ],
                  resid[lr_base$receptor[i], ])
  row <- tibble(
    ligand = lr_base$ligand[i], receptor = lr_base$receptor[i],
    pooled_rho = pooled[1], pooled_pvalue = pooled[2],
    residual_rho = adjusted[1], residual_pvalue = adjusted[2])
  for (lev in unique(condition))
    row[[paste0("within_", make.names(lev), "_rho")]] <-
      rho(l[condition == lev], r[condition == lev])[1]
  row
})) |>
  mutate(pooled_padj = p.adjust(pooled_pvalue, "BH"),
         residual_padj = p.adjust(residual_pvalue, "BH"))

de_list <- readRDS(obj_path("de_list.rds"))
for (name in names(de_list)) {
  sig <- de_list[[name]]$gene[de_list[[name]]$significant]
  pair_stats[[paste0("associated_", name)]] <-
    pair_stats$ligand %in% sig | pair_stats$receptor %in% sig
}
assoc <- grep("^associated_", names(pair_stats), value = TRUE)
vs_control <- Filter(
  function(x) identical(x$denominator, cfg$design$reference_level),
  cfg$comparisons)
tgf_cmp <- Filter(function(x) identical(x$numerator, "TGF"), vs_control)
bleo_cmp <- Filter(function(x) identical(x$numerator, "Bleo"), vs_control)
tgf_col <- if (length(tgf_cmp)) paste0("associated_", tgf_cmp[[1]]$name) else ""
bleo_col <- if (length(bleo_cmp)) paste0("associated_", bleo_cmp[[1]]$name) else ""
in_t <- if (tgf_col %in% assoc) pair_stats[[tgf_col]] else FALSE
in_b <- if (bleo_col %in% assoc) pair_stats[[bleo_col]] else FALSE
pair_stats$associated_in <- case_when(
  in_t & in_b ~ "Both treatments",
  in_b ~ display_condition("Bleo"),
  in_t ~ display_condition("TGF"),
  TRUE ~ "None")

expand_candidates <- function(membership) {
  lig <- membership |>
    select(ligand = gene, ligand_ct = cell_type,
           ligand_reference_expression = reference_expression)
  rec <- membership |>
    select(receptor = gene, receptor_ct = cell_type,
           receptor_reference_expression = reference_expression)
  lr_base |>
    inner_join(lig, by = "ligand", relationship = "many-to-many") |>
    inner_join(rec, by = "receptor", relationship = "many-to-many") |>
    filter(ligand_ct != receptor_ct) |>
    mutate(reference_score = sqrt(ligand_reference_expression *
                                    receptor_reference_expression)) |>
    inner_join(pair_stats, by = c("ligand", "receptor")) |>
    arrange(ligand_ct, receptor_ct, desc(reference_score), ligand, receptor)
}

# --- primary: five GSE136831 states ----------------------------------------
profiles <- load_GSE136831_profiles()
membership <- profiles$membership |>
  filter(member) |>
  transmute(gene, cell_type, reference_expression = mean_cpm)
pairs <- expand_candidates(membership)
if (!nrow(pairs))
  stop("No GSE136831 five-state candidate L-R pairs passed the thresholds.",
       call. = FALSE)

state_levels <- profiles$states
pair_counts <- tidyr::expand_grid(ligand_ct = state_levels,
                                  receptor_ct = state_levels) |>
  filter(ligand_ct != receptor_ct) |>
  left_join(count(pairs, ligand_ct, receptor_ct, name = "n_pairs"),
            by = c("ligand_ct", "receptor_ct")) |>
  mutate(n_pairs = coalesce(n_pairs, 0L))
message("Primary GSE136831 five-state candidate pairs: ", nrow(pairs))
print(pair_counts)
write_csv(pairs, tbl_path("lr_candidates.csv"))
saveRDS(pairs, obj_path("lr_candidates.rds"))

# --- treatment-associated candidates ---------------------------------------
# The chord panels come from Stage 08; this is the only panel that splits the
# five-state candidates by which contrast made the ligand or receptor
# differentially expressed, so it stays here. Labels come from the config
# display map, the colours are the original panel's.
treatment_levels <- c(display_condition("TGF"), display_condition("Bleo"),
                      "Both treatments")
treatment <- pairs |>
  filter(associated_in != "None") |>
  slice_max(reference_score, n = 30) |>
  mutate(pair = paste(ligand, "→", receptor),
         axis = paste(ligand_ct, "→", receptor_ct),
         associated_in = factor(associated_in, treatment_levels))
save_plot(
  ggplot(treatment, aes(axis, reorder(pair, reference_score),
                        colour = associated_in, size = reference_score)) +
    geom_point() +
    scale_colour_manual(
      values = setNames(c("#4DBBD5", "#E64B35", "#3C5488"), treatment_levels),
      drop = FALSE) +
    labs(x = "sender → receiver", y = NULL, colour = "associated in",
         size = "reference score",
         title = "Inferenced cell-cell communication") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1)),
  "lr_treatment_associated.svg", 8, 6)

message("07 done -> ", OUT)
