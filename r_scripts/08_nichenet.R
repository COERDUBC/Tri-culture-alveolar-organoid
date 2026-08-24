#!/usr/bin/env Rscript
# 08 — NicheNet ligand-activity prioritisation.
#
# Different question from Stage 07: instead of asking "are both partners
# expressed?", NicheNet asks "do a ligand's known downstream targets recover the
# genes that actually changed under treatment?" That makes it a functional
# ranking of candidate ligands, using prior regulatory networks.
#
# Needs the prior model files (config: nichenet$ligand_target_matrix, lr_network).
# Skips cleanly when nichenetr or the files are unavailable.

source("r_scripts/common.R")
suppressPackageStartupMessages({ library(edgeR) })

ltm_path <- cfg$nichenet$ligand_target_matrix %||% ""
lr_path  <- cfg$nichenet$lr_network %||% ""
ok <- function(p) nzchar(p) && !grepl("^TO_FILL$", p) && file.exists(p)
if (!requireNamespace("nichenetr", quietly = TRUE) || !ok(ltm_path) || !ok(lr_path)) {
  message("08 skipped — needs nichenetr plus nichenet$ligand_target_matrix and ",
          "nichenet$lr_network in config.yml.")
  quit(save = "no")
}
if (!requireNamespace("circlize", quietly = TRUE))
  stop("Stage 08 requires `circlize`; run Rscript setup_packages.R first.",
       call. = FALSE)
m <- readRDS(obj_path("model.rds"))
de_list <- readRDS(obj_path("de_list.rds"))
top_n <- cfg$nichenet$top_n_ligands %||% 20
top_n_chord <- cfg$nichenet$top_n_chord_pairs %||% 40
lr_primary_path <- obj_path("lr_candidates.rds")
if (!file.exists(lr_primary_path))
  stop("Stage 08 requires the GSE136831 five-state Stage 07 result: ",
       lr_primary_path, call. = FALSE)
lr_primary <- readRDS(lr_primary_path)
profiles <- load_GSE136831_profiles()
state_membership <- profiles$membership |>
  filter(member) |>
  select(gene, cell_type)
ligand_annotations <- state_membership |>
  group_by(gene) |>
  summarise(
    ligand_states = paste(sort(unique(cell_type)), collapse = ";"),
    n_ligand_states = n_distinct(cell_type),
    .groups = "drop")

ltm <- readRDS(ltm_path)                       # ligand x target regulatory potential
lr_network <- readRDS(lr_path)                 # from / to
message("Prior model: ", nrow(ltm), " targets x ", ncol(ltm), " ligands")

cpm_lin <- cpm(to_symbol(m$dge$counts), log = FALSE)
expressed <- rownames(cpm_lin)[rowMeans(cpm_lin >= 1) >= 0.5]

results <- list()
for (name in names(de_list)) {
  d <- de_list[[name]]
  # Gene set of interest = what changed; background = what was testable.
  background <- intersect(expressed, rownames(ltm))
  geneset <- intersect(unique(d$gene[d$significant]), background)
  if (length(geneset) < 10) { message(name, ": too few DE genes — skipped"); next }

  # Stage 07 has already required expression and a receptor in a different
  # GSE136831 state. Also require an expressed receptor in the NicheNet prior.
  prior_ligands <- lr_network |>
    filter(from %in% expressed, to %in% expressed) |>
    pull(from) |>
    unique()
  ligands <- unique(lr_primary$ligand) |>
    intersect(prior_ligands) |>
    intersect(colnames(ltm))
  if (!length(ligands)) { message(name, ": no candidate ligand — skipped"); next }

  act <- nichenetr::predict_ligand_activities(
    geneset = geneset, background_expressed_genes = background,
    ligand_target_matrix = ltm, potential_ligands = ligands) |>
    arrange(desc(aupr_corrected)) |>
    mutate(rank = row_number()) |>
    left_join(ligand_annotations, by = c("test_ligand" = "gene")) |>
    mutate(
      ligand_cell_type = if_else(
        n_ligand_states == 1L, ligand_states, "Multiple states"),
      comparison = name, .before = 1)
  write_csv(act, tbl_path(paste0("nichenet_ligand_activity_", name, ".csv")))
  results[[name]] <- act
  message(name, ": ", nrow(act), " ligands ranked; top = ",
          paste(head(act$test_ligand, 5), collapse = ", "))

  p <- act |> slice_head(n = top_n) |>
    ggplot(aes(aupr_corrected, reorder(test_ligand, aupr_corrected),
               fill = ligand_cell_type)) +
    geom_col() + scale_fill_manual(values = PAL, na.value = "grey75",
                                   name = "ligand cell type") +
    labs(title = paste("NicheNet ligand activity:", display_comparison(name)),
         x = "AUPR (corrected)", y = NULL) + theme_bw(base_size = 9)
  save_plot(p, paste0("nichenet_lr_activity_", name, ".svg"), 7, 6)
}

if (length(results)) saveRDS(results, obj_path("nichenet.rds"))

# Candidate L-R pairs whose ligand is among the top NicheNet-ranked ligands.
# The highest reference-score subset is shown as a directional chord; the
# complete intersected table is retained for audit and interpretation.
if (length(results)) {
  lr <- lr_primary
  present <- unique(c(lr$ligand_ct, lr$receptor_ct))
  preferred <- as.character(unlist(
    cfg$deconvolution$GSE136831_atlas$lineage_restricted_cell_types))
  ct_levels <- c(intersect(preferred, present),
                 sort(setdiff(present, preferred)))
  ct_col <- setNames(PAL[seq_along(ct_levels)], ct_levels)

  for (name in names(results)) {
    active <- results[[name]] |>
      slice_head(n = top_n) |>
      transmute(
        ligand = test_ligand,
        ligand_rank = rank,
        ligand_aupr_corrected = aupr_corrected)
    sub <- lr |>
      inner_join(active, by = "ligand")
    if (!nrow(sub)) {
      message(name, ": no candidate pair has a top-", top_n, " active ligand.")
      next
    }
    write_csv(sub, tbl_path(paste0("nichenet_active_lr_", name, ".csv")))
    display_pairs <- sub |>
      filter(is.finite(reference_score), reference_score > 0) |>
      slice_max(reference_score, n = top_n_chord, with_ties = FALSE)
    if (!nrow(display_pairs)) {
      message(name, ": all displayed pair widths were zero.")
      next
    }
    df <- display_pairs |>
      transmute(
        from = paste(ligand_ct, ligand, sep = ":"),
        to = paste(receptor_ct, receptor, sep = ":"),
        value = log2(reference_score + 1))

    ligand_nodes <- paste(
      display_pairs$ligand_ct, display_pairs$ligand, sep = ":")
    receptor_nodes <- paste(
      display_pairs$receptor_ct, display_pairs$receptor, sep = ":")
    gene_ct <- c(
      setNames(display_pairs$ligand_ct, ligand_nodes),
      setNames(display_pairs$receptor_ct, receptor_nodes))
    gene_ct <- gene_ct[!duplicated(names(gene_ct))]
    save_lr_flow(
      df, gene_ct, ct_col,
      paste0("nichenet_lr_chord_", name, ".svg"),
      paste0("Active-ligand L-R: ", display_comparison(name)),
      show_cell_type_prefix = FALSE,
      inner_score = display_pairs$reference_score)

    message(name, ": ", nrow(sub), " candidate pairs from ",
            n_distinct(sub$ligand), " active ligands; ", nrow(display_pairs),
            " highest reference-score pairs shown in the gene chord")
  }
}

message("08 done -> ", OUT)
