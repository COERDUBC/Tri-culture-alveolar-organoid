#!/usr/bin/env Rscript
# 03 — Ranked gene-set enrichment. Three layers:
#   (a) curated IPF/fibrosis panels  = the prespecified, interpretable test
#   (b) config-selected MSigDB collections = broad pathway safety net
#   (c) published Ryu 2024 aberrant-basaloid IPF-severity programme
# Genes are ranked by the moderated t-statistic, so ranking matches the DE test.

source("r_scripts/common.R")
suppressPackageStartupMessages({ library(fgsea); library(msigdbr) })

de_list <- readRDS(obj_path("de_list.rds"))
configured_comparisons <- vapply(
  cfg$comparisons, function(x) as.character(x$name), character(1))
if (!setequal(names(de_list), configured_comparisons))
  stop("de_list.rds comparisons do not match config.yml. Re-run Stage 02 ",
       "for this run before Stage 03.", call. = FALSE)

# Curated panels. TGF-β is expected to drive Myofibroblast_ECM; Bleomycin the
# Aberrant_basaloid / Senescence axis (Adams 2020, Habermann 2020, Tsukui 2020).
panels <- list(
  AT1 = c("AGER","AQP5","HOPX","PDPN","CAV1","RTKN2","CLIC5","EMP2"),
  AT2 = c("SFTPC","SFTPA1","SFTPB","SFTPD","LAMP3","ABCA3","SLC34A2","NAPSA"),
  Epithelial_polarity = c("CDH1","OCLN","EPCAM","TJP1","CLDN3","CLDN4","CLDN7","KRT8","KRT18"),
  Fibroblast = c("COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","PDGFRB","FN1","SPARC","ELN"),
  Myofibroblast_ECM = c("ACTA2","COL1A1","COL1A2","COL3A1","CTHRC1","FN1","POSTN",
                        "SERPINE1","TAGLN","FBN1","TIMP1","LOX"),
  Aberrant_basaloid = c("KRT17","KRT8","CDKN1A","CDKN2A","ITGB6","GDF15","EPHB2","SFN","MMP7"),
  Senescence_SASP = c("CDKN1A","CDKN2A","TP53","GDF15","IGFBP3","IGFBP7","SERPINE1",
                      "IL1A","IL1B","IL6","CXCL8","CCL2","MMP1","MMP3"),
  EMT = c("VIM","FN1","CDH2","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","ITGA5","MMP2","MMP9"))

# MSigDB renamed the KEGG subcollections (CP:KEGG -> CP:KEGG_MEDICUS /
# CP:KEGG_LEGACY), so try the current names in order.
kegg_sets <- function() {
  for (sub in c("CP:KEGG_MEDICUS", "CP:KEGG_LEGACY", "CP:KEGG")) {
    s <- tryCatch(msig("C2", sub), error = function(e) NULL)
    if (length(s)) { message("  KEGG subcollection: ", sub, " (", length(s), " sets)"); return(s) }
  }
  list()
}

# --- configured GSEA collections -------------------------------------------
# User-facing GSEA collection registry. The config controls computation and
# output together; curated panels and the separate Ryu module are unaffected.
# The shape of this config block (ids, supported sources and collections,
# top_terms_per_plot) is set once in config.yml.
gene_set_cfg <- cfg$gsea$gene_sets %||% list()
if (!length(gene_set_cfg))
  stop("gsea$gene_sets must contain at least one configured collection.",
       call. = FALSE)
gene_set_specs <- bind_rows(lapply(gene_set_cfg, function(x) tibble(
  id = as.character(x$id),
  source = as.character(x$source),
  collection = as.character(x$collection))))

top_terms_per_plot <- as.integer(cfg$gsea$top_terms_per_plot %||% 15L)

# --- load collections ------------------------------------------------------
load_collection <- function(collection) switch(
  collection,
  Hallmark = msig("H"),
  GO_BP = msig("C5", "GO:BP"),
  Reactome = msig("C2", "CP:REACTOME"),
  KEGG = kegg_sets())
message("Loading configured MSigDB collections: ",
        paste(gene_set_specs$collection, collapse = ", "))
collections <- setNames(
  lapply(gene_set_specs$collection, load_collection),
  gene_set_specs$id)
missing_collections <- names(collections)[lengths(collections) == 0L]
if (length(missing_collections))
  stop("Configured MSigDB collection(s) unavailable: ",
       paste(paste0(missing_collections, " (",
                    gene_set_specs$collection[
                      match(missing_collections, gene_set_specs$id)], ")"),
             collapse = ", "), call. = FALSE)
message("Configured gene sets: ",
        paste(sprintf("%s=%s (%d)", gene_set_specs$id,
                      gene_set_specs$collection, lengths(collections)),
              collapse = ", "))

# --- collection dotplot ----------------------------------------------------
# Overview dotplot for every configured collection. `id` controls the stable
# filename while `collection` records the MSigDB family used for computation.
# Called once per contrast from the loop below, which supplies `name` and the
# contrast's `gsea_res`.
collection_strip <- c(
  Hallmark = "^HALLMARK_",
  GO_BP = "^(GOBP_|GO_)",
  Reactome = "^REACTOME_",
  KEGG = "^KEGG_(MEDICUS_|LEGACY_)?")
dotplot_collection <- function(id, name, gsea_res) {
  meta <- gene_set_specs[match(id, gene_set_specs$id), ]
  strip <- unname(collection_strip[[meta$collection]])
  d <- gsea_res |> filter(gene_set_id == id) |>
    slice_min(padj, n = top_terms_per_plot, with_ties = FALSE) |>
    mutate(label = gsub("_", " ", sub(strip, "", pathway)),
           plot_padj = pmax(padj, .Machine$double.xmin))
  if (!nrow(d))
    stop("No GSEA result available for configured collection ", id,
         " (", meta$collection, ") in contrast ", name, ".",
         call. = FALSE)
  p <- ggplot(d, aes(NES, reorder(label, NES), size = size, colour = plot_padj)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point() +
    scale_colour_gradient(low = "#E64B35", high = "#4DBBD5", trans = "log10") +
    scale_size_continuous(range = c(2, 7)) +
    labs(title = paste0(gsub("_", " ", meta$collection), " GSEA: ",
                        display_comparison(name)),
         y = NULL, x = "NES",
         size = "set size", colour = "adjusted P") +
    theme_bw(base_size = 9)
  save_plot(
    p, paste0("gsea_", id, "_dotplot_", name, ".svg"), 7.5, 5)
}

# --- fgsea helpers ---------------------------------------------------------
run <- function(pathways, ranks, min_size) {
  set.seed(SEED)
  fgsea(pathways, ranks, minSize = min_size, maxSize = 500, eps = 0) |>
    as_tibble() |> arrange(padj)
}
# One fgsea call per collection, so fgsea's BH is computed within that collection.
run_by_collection <- function(ranks, min_size) {
  bind_rows(lapply(names(collections), function(id) {
    meta <- gene_set_specs[match(id, gene_set_specs$id), ]
    run(collections[[id]], ranks, min_size) |>
      mutate(
        gene_set_id = id,
        source = meta$source,
        collection = meta$collection,
        .before = 1)
  })) |>
    arrange(gene_set_id, padj)
}

# --- Ryu 2024 module gene set ----------------------------------------------
# Kept separate from the hand-curated marker panels above: those are BH-adjusted
# within their eight-panel family, and folding this set in would shift every one
# of their adjusted P values. Listed in config.yml rather than read from a GMT,
# like the other gene lists in this project.
ryu_cfg <- cfg$published_signatures$ryu_2024_module
ryu_name <- as.character(ryu_cfg$name %||% "")
ryu_source <- unique(trimws(as.character(unlist(ryu_cfg$genes))))
if (!nzchar(ryu_name) || !length(ryu_source) || any(!nzchar(ryu_source)))
  stop("published_signatures$ryu_2024_module needs a name and a non-empty ",
       "genes list.", call. = FALSE)

panel_summaries <- list()
msigdb_summaries <- list()
ryu_summaries <- list()
ryu_leading_edges <- list()
coverage_summaries <- list()
uncovered_gene_tables <- list()
leading_edge_maps <- list()
curated_panel_genes <- unique(unname(unlist(panels, use.names = FALSE)))
coverage_top_n <- as.integer(cfg$panel_coverage$top_n_de_genes %||% 50L)
coverage_msigdb_padj <- cfg$panel_coverage$msigdb_padj %||% cfg$thresholds$padj

# --- per-contrast GSEA -----------------------------------------------------
for (name in names(de_list)) {
  d <- de_list[[name]] |> filter(!is.na(gene), gene != "", is.finite(t))
  # One value per symbol: keep the most extreme moderated t.
  rank_table <- d |>
    arrange(gene, desc(abs(t))) |>
    distinct(gene, .keep_all = TRUE)
  ranks <- setNames(rank_table$t, rank_table$gene)
  message("GSEA ", name, ": ", length(ranks), " ranked symbols")

  panel_res <- run(panels, ranks, 3)
  gsea_res  <- run_by_collection(ranks, 10)
  panel_summaries[[name]] <- panel_res |>
    mutate(comparison = name, .before = 1)
  msigdb_summaries[[name]] <- gsea_res |>
    mutate(comparison = name, .before = 1)

  # Coverage audit of the strict top-DE genes. The comparison set is exactly
  # the eight hand-curated panels above; the separately tested Ryu programme is
  # deliberately not counted as panel coverage. Pathway annotations require a
  # gene to be in the leading edge of a significant MSigDB term, rather than
  # merely being listed somewhere in a large database gene set.
  top_de <- rank_table |>
    filter(significant) |>
    arrange(desc(abs(t)), padj, gene) |>
    slice_head(n = coverage_top_n) |>
    mutate(
      comparison = name,
      comparison_label = display_comparison(name),
      top_de_rank = row_number(),
      covered_by_curated_panel = gene %in% curated_panel_genes,
      .before = 1)
  uncovered <- top_de |>
    filter(!covered_by_curated_panel) |>
    select(comparison, comparison_label, top_de_rank, gene, gene_id,
           log2FC, t, pvalue, padj, direction)

  significant_msigdb <- gsea_res |>
    filter(!is.na(padj), padj < coverage_msigdb_padj)
  leading_edge_map <- tibble(
    comparison = character(), gene = character(), gene_set_id = character(),
    source = character(), collection = character(),
    pathway = character(), pathway_NES = double(), pathway_padj = double(),
    pathway_size = integer())
  if (nrow(significant_msigdb)) {
    leading_edge_map <- bind_rows(lapply(seq_len(nrow(significant_msigdb)),
      function(i) {
        genes <- unique(as.character(significant_msigdb$leadingEdge[[i]]))
        tibble(
          comparison = name, gene = genes,
          gene_set_id = significant_msigdb$gene_set_id[[i]],
          source = significant_msigdb$source[[i]],
          collection = significant_msigdb$collection[[i]],
          pathway = significant_msigdb$pathway[[i]],
          pathway_NES = significant_msigdb$NES[[i]],
          pathway_padj = significant_msigdb$padj[[i]],
          pathway_size = significant_msigdb$size[[i]])
      }))
  }
  leading_edge_maps[[name]] <- leading_edge_map
  uncovered_pathways <- uncovered |>
    inner_join(leading_edge_map, by = c("comparison", "gene")) |>
    arrange(top_de_rank, pathway_padj, collection, pathway)
  pathway_summary <- uncovered_pathways |>
    group_by(comparison, gene) |>
    summarise(
      n_msigdb_leading_edge_pathways = n(),
      msigdb_collections = paste(sort(unique(collection)), collapse = ";"),
      msigdb_leading_edge_pathways = paste(
        paste(collection, pathway, sep = ":"), collapse = ";"),
      .groups = "drop")
  uncovered <- uncovered |>
    left_join(pathway_summary, by = c("comparison", "gene")) |>
    mutate(
      n_msigdb_leading_edge_pathways = coalesce(
        n_msigdb_leading_edge_pathways, 0L),
      msigdb_collections = coalesce(msigdb_collections, "None"),
      msigdb_leading_edge_pathways = coalesce(
        msigdb_leading_edge_pathways, "None"))

  coverage_summaries[[name]] <- tibble(
    comparison = name,
    comparison_label = display_comparison(name),
    strict_de_genes = sum(rank_table$significant),
    top_de_genes_audited = nrow(top_de),
    covered_by_curated_panels = sum(top_de$covered_by_curated_panel),
    uncovered_from_curated_panels = nrow(uncovered),
    uncovered_with_significant_msigdb_leading_edge =
      sum(uncovered$n_msigdb_leading_edge_pathways > 0),
    curated_panel_coverage_fraction = if (nrow(top_de))
      mean(top_de$covered_by_curated_panel) else NA_real_)
  uncovered_gene_tables[[name]] <- uncovered

  ryu_genes <- intersect(ryu_source, names(ranks))
  ryu_res <- run(setNames(list(ryu_genes), ryu_name), ranks, ryu_cfg$min_size)
  ryu_summaries[[name]] <- ryu_res |>
    mutate(comparison = name, source_set_size = length(ryu_source),
           matched_set_size = length(ryu_genes),
           source_table = as.character(ryu_cfg$source_table),
           doi = as.character(ryu_cfg$doi), .before = 1) |>
    select(-padj)

  # Preserve one auditable, combined leading-edge table across contrasts.
  ryu_gene_table <- tibble(gene = unique(ryu_res$leadingEdge[[1]])) |>
    left_join(rank_table |>
                select(any_of(c("gene", "t", "log2FC", "pvalue", "padj",
                                "significant", "direction"))),
              by = "gene") |>
    mutate(
      comparison = name,
      programme = ryu_name,
      programme_NES = ryu_res$NES[[1]],
      programme_pvalue = ryu_res$pval[[1]],
      enrichment_direction = if_else(ryu_res$NES[[1]] > 0,
                                     "up_with_treatment", "down_with_treatment"),
      .before = 1) |>
    arrange(desc(abs(t)), gene)
  ryu_leading_edges[[name]] <- ryu_gene_table

  # Curated panels as a bar chart — this is the headline confirmatory result.
  p <- panel_res |> mutate(pathway = reorder(pathway, NES)) |>
    ggplot(aes(NES, pathway, fill = padj < cfg$thresholds$padj)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "grey75"),
                      name = paste0("padj < ", cfg$thresholds$padj)) +
    labs(title = paste("Curated panels:", display_comparison(name)), y = NULL) +
    theme_bw()
  save_plot(p, paste0("panels_", name, ".svg"), 7, 4)

  for (id in gene_set_specs$id) dotplot_collection(id, name, gsea_res)

}

# --- combined tables -------------------------------------------------------
# A configured collection is an output contract: every configured DE contrast
# must have its collection-specific dot plot before Stage 03 can finish.
expected_gsea_figures <- as.vector(outer(
  gene_set_specs$id, names(de_list),
  function(id, comparison)
    fig_path(paste0("gsea_", id, "_dotplot_", comparison, ".svg"))))
missing_gsea_figures <- expected_gsea_figures[!file.exists(expected_gsea_figures)]
if (length(missing_gsea_figures))
  stop("Stage 03 did not create every configured GSEA dot plot: ",
       paste(basename(missing_gsea_figures), collapse = ", "), call. = FALSE)

write_csv(flatten_le(bind_rows(panel_summaries) |>
                       arrange(comparison, padj)),
          tbl_path("gsea_panels.csv"))
write_csv(flatten_le(bind_rows(msigdb_summaries) |>
                       arrange(comparison, gene_set_id, padj)),
          tbl_path("gsea_msigdb.csv"))
write_csv(bind_rows(coverage_summaries),
          tbl_path("curated_panel_coverage_summary.csv"))
write_csv(bind_rows(uncovered_gene_tables),
          tbl_path("curated_panel_uncovered_top_de_genes.csv"))

# --- leading-edge overlap between enriched sets -----------------------------
# MSigDB sets are not mutually exclusive, so several co-enriched terms can be
# driven by the same genes and must not be read as independent findings. The
# sets are the ones already filtered at `coverage_msigdb_padj` for the coverage
# audit, so both analyses quote the same enriched terms. Both tables are
# computed within a comparison, collection and enrichment direction:
# sharing between an up- and a down-enriched set is a different question and
# would only dilute the redundancy the reader is actually asking about.
leading_edge <- bind_rows(leading_edge_maps) |>
  mutate(direction = if_else(pathway_NES > 0, "up", "down"))

# `n_pathways_sharing_gene == 1` gives each set's leading-edge genes that no
# other enriched set in the same direction claims — the specificity number to
# quote alongside an enrichment.
write_csv(
  leading_edge |>
    group_by(comparison, collection, direction, gene) |>
    mutate(n_pathways_sharing_gene = n()) |>
    ungroup() |>
    arrange(comparison, collection, pathway_padj, pathway, gene),
  tbl_path("gsea_leading_edge_membership.csv"))

# Pairs with no shared leading-edge gene are the uninformative majority and are
# dropped, which is what keeps the Reactome contrasts to a readable table.
leading_edge_overlap <- leading_edge |>
  group_by(comparison, collection, direction) |>
  group_modify(function(d, key) {
    sets <- split(d$gene, d$pathway)
    if (length(sets) < 2) return(tibble())
    nes <- deframe(distinct(d, pathway, pathway_NES))
    combos <- combn(names(sets), 2)
    tibble(
      pathway_a = combos[1, ], pathway_b = combos[2, ],
      NES_a = unname(nes[combos[1, ]]), NES_b = unname(nes[combos[2, ]]),
      n_leading_edge_a = unname(lengths(sets[combos[1, ]])),
      n_leading_edge_b = unname(lengths(sets[combos[2, ]])),
      shared_genes = vapply(seq_len(ncol(combos)), function(i)
        paste(sort(intersect(sets[[combos[1, i]]], sets[[combos[2, i]]])),
              collapse = ";"), character(1))) |>
      filter(nzchar(shared_genes)) |>
      mutate(
        n_shared = lengths(strsplit(shared_genes, ";")),
        jaccard = n_shared /
          (n_leading_edge_a + n_leading_edge_b - n_shared))
  }) |>
  ungroup() |>
  arrange(comparison, collection, desc(jaccard))
write_csv(leading_edge_overlap, tbl_path("gsea_leading_edge_overlap.csv"))
message("Leading-edge overlap: ", nrow(leading_edge_overlap),
        " enriched set pairs share at least one leading-edge gene")

ryu_summary <- bind_rows(ryu_summaries) |>
  arrange(comparison) |>
  mutate(comparison_label = factor(display_comparison(comparison),
                                    levels = rev(display_comparison(c("TGF_vs_Control", "Bleo_vs_Control", "TGF_vs_Bleo")))))
write_csv(flatten_le(ryu_summary), tbl_path("gsea_ryu2024.csv"))
ryu_leading_edge <- bind_rows(ryu_leading_edges)
write_csv(ryu_leading_edge, tbl_path("ryu2024_leading_edge.csv"))

# --- Ryu A/B composite figure ----------------------------------------------
ryu_summary_plot <-
  ggplot(ryu_summary,
         aes(NES, comparison_label, fill = NES > 0)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_col(width = 0.65) +
  geom_text(aes(label = paste0("P=", signif(pval, 2))),
            hjust = ifelse(ryu_summary$NES > 0, -0.08, 1.08), size = 3) +
  scale_fill_manual(values = c(`TRUE` = "#E64B35", `FALSE` = "#4DBBD5"),
                    guide = "none") +
  scale_x_continuous(expand = expansion(mult = 0.12)) +
  # No title or subtitle: the gene set is identified in the manuscript figure
  # legend, not on the panel.
  labs(
    tag = "A",
    x = "NES", y = NULL) +
  theme_bw(base_size = 9) +
  theme(plot.tag = element_text(face = "bold", size = 13))

# Compact gene-level companion to the programme NES summary. 
ryu_heatmap_members <- bind_rows(lapply(names(de_list), function(name) {
  de_list[[name]] |>
    filter(gene %in% ryu_source, is.finite(t)) |>
    arrange(gene, desc(abs(t))) |>
    distinct(gene, .keep_all = TRUE) |>
    transmute(comparison = name, gene, t)
}))
ryu_heatmap_genes <- unique(ryu_heatmap_members$gene)

if (length(ryu_heatmap_genes)) {
  ryu_heatmap_effects <- bind_rows(lapply(names(de_list), function(name) {
    de_list[[name]] |>
      filter(!is.na(gene), gene %in% ryu_heatmap_genes, is.finite(t)) |>
      arrange(gene, desc(abs(t))) |>
      distinct(gene, .keep_all = TRUE) |>
      transmute(comparison = name, gene, log2FC, t, padj)
  }))
  ryu_heatmap_grid <- tidyr::expand_grid(
    comparison = configured_comparisons,
    gene = ryu_heatmap_genes) |>
    left_join(ryu_heatmap_effects, by = c("comparison", "gene")) |>
    left_join(
      ryu_leading_edge |>
        distinct(comparison, gene) |>
        mutate(leading_edge = TRUE),
      by = c("comparison", "gene")) |>
    mutate(
      leading_edge = coalesce(leading_edge, FALSE),
      comparison_label = display_comparison(comparison))
  if (any(!is.finite(ryu_heatmap_grid$log2FC)))
    stop("Ryu leading-edge heatmap could not recover a finite log2FC for ",
         "every selected gene and contrast.", call. = FALSE)

  gene_order <- ryu_heatmap_members |>
    mutate(comparison_order = match(comparison, configured_comparisons),
           abs_t = abs(t)) |>
    group_by(gene) |>
    slice_max(abs_t, n = 1, with_ties = FALSE) |>
    ungroup() |>
    arrange(comparison_order, desc(abs_t), gene) |>
    pull(gene)
  ryu_heatmap_grid <- ryu_heatmap_grid |>
    mutate(
      gene = factor(gene, levels = rev(unique(gene_order))),
      comparison_label = factor(
        comparison_label,
        levels = display_comparison(configured_comparisons))) |>
    arrange(gene, comparison_label)
  write_csv(
    ryu_heatmap_grid |>
      mutate(gene = as.character(gene),
             comparison_label = as.character(comparison_label)),
    tbl_path("ryu2024_leading_edge_heatmap.csv"))

  heatmap_height <- max(4.5, min(10, 2 + 0.22 * length(gene_order)))
  ryu_heatmap_plot <-
    ggplot(ryu_heatmap_grid,
           aes(comparison_label, gene, fill = log2FC)) +
      geom_tile(colour = "white", linewidth = 0.3) +
      geom_point(
        data = filter(ryu_heatmap_grid, leading_edge),
        shape = 21, fill = NA, colour = "black", stroke = 0.5, size = 2.3) +
      diverging_fill(ryu_heatmap_grid$log2FC) +
      # No title, subtitle or caption: the gene set and the circle marker are
      # both explained in the manuscript figure legend, not on the panel.
      labs(
        tag = "B",
        x = NULL, y = NULL, fill = "limma log2FC") +
      theme_bw(base_size = 9) +
      theme(
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 30, hjust = 1),
        axis.text.y = element_text(face = "italic"),
        plot.tag = element_text(face = "bold", size = 13))

  ryu_summary_plot <- prepare_plot_for_export(ryu_summary_plot)
  ryu_heatmap_plot <- prepare_plot_for_export(ryu_heatmap_plot)
  save_figure(
    "gsea_ryu2024_programme_and_leading_edge.svg", 8.5,
    3.4 + heatmap_height,
    function() {
      grid::grid.newpage()
      ryu_layout <- grid::grid.layout(
        nrow = 2, ncol = 1,
        heights = grid::unit(c(3.4, heatmap_height), "null"))
      grid::pushViewport(grid::viewport(layout = ryu_layout))
      print(ryu_summary_plot,
            vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
      print(ryu_heatmap_plot,
            vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
      grid::popViewport()
    })
}
if (!file.exists(fig_path("gsea_ryu2024_programme_and_leading_edge.svg")) ||
    !file.exists(tbl_path("ryu2024_leading_edge_heatmap.csv")))
  stop("Stage 03 did not create the Ryu A/B figure and plotting table.",
       call. = FALSE)

message("03 done -> ", OUT)
