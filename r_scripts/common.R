# Shared helpers used by multiple analysis stages.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(ggplot2)
})

# config.yml is the committed template and carries <ProjectPath> placeholders.
# config_local.yml holds this machine's real paths and is git-ignored, so it
# wins when present; CONFIG overrides both.
CONFIG_FILE <- Sys.getenv("CONFIG", "")
if (!nzchar(CONFIG_FILE))
  CONFIG_FILE <- if (file.exists("config_local.yml")) "config_local.yml" else "config.yml"
if (!file.exists(CONFIG_FILE))
  stop("Configuration file not found: ", CONFIG_FILE, call. = FALSE)
message("Config: ", CONFIG_FILE)
cfg <- yaml::read_yaml(CONFIG_FILE)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Analysis seed, resolved once so the default cannot drift between stages.
# 02_de.R keeps its own literal jitter seed: that one positions points in a
# figure and should not move when the analysis seed changes.
SEED <- cfg$reproducibility$random_seed %||% 20260715

# "Missing: a, b; extra: c" — the shared tail of every set-mismatch stop().
set_diff_text <- function(expected, observed)
  paste0("Missing: ", paste(setdiff(expected, observed), collapse = ", "),
         "; extra: ", paste(setdiff(observed, expected), collapse = ", "))

# One output folder per pipeline run. Stage 01 requests a new folder; later
# stages inherit it through <results_dir>/CURRENT_RUN. Setting OUT
# explicitly always wins, which permits intentional resume/re-runs.
out_dir <- function(start_new = isTRUE(getOption("fibres.start_new_run", FALSE))) {
  root <- cfg$output$results_dir %||% "results"
  pointer <- file.path(root, "CURRENT_RUN")
  run <- Sys.getenv("OUT", "")
  if (!nzchar(run) && start_new) {
    base <- format(Sys.time(), "run_%Y-%m-%d_%H-%M-%S")
    run <- base
    suffix <- 1L
    while (dir.exists(file.path(root, run))) {
      run <- paste0(base, "_", suffix)
      suffix <- suffix + 1L
    }
  }
  if (!nzchar(run) && file.exists(pointer)) run <- trimws(readLines(pointer, n = 1))
  if (!nzchar(run))
    stop("No active analysis run. Run Stage 01 first or set OUT to an existing ",
         "run directory name.", call. = FALSE)
  if (!grepl("^run_[A-Za-z0-9._-]+$", run))
    stop("OUT/CURRENT_RUN must be a safe run_* directory name: ", run,
         call. = FALSE)
  if (!start_new && !dir.exists(file.path(root, run)))
    stop("Active analysis run does not exist: ", file.path(root, run),
         ". Run Stage 01 first or set OUT correctly.", call. = FALSE)
  dir.create(file.path(root, run), recursive = TRUE, showWarnings = FALSE)
  writeLines(run, pointer)
  file.path(root, run)
}

OUT <- out_dir()
for (d in c("tables", "figures", "objects"))
  dir.create(file.path(OUT, d), recursive = TRUE, showWarnings = FALSE)

# Refresh the run's configuration after both pre- and post-web stages.
file.copy(CONFIG_FILE, file.path(OUT, "config.yml"), overwrite = TRUE)
tracked_packages <- c(
  "tximport", "edgeR", "limma", "fgsea", "msigdbr", "circlize", "nichenetr")
readr::write_csv(
  tibble::tibble(
    package = tracked_packages,
    # packageVersion() reads DESCRIPTION; requireNamespace() would LOAD each package,
    # which pulled nichenetr (and its pROC import) into every stage and buried the
    # real warnings under its startup messages.
    version = vapply(tracked_packages, function(pkg)
      tryCatch(as.character(utils::packageVersion(pkg)),
               error = function(e) NA_character_), character(1))),
  file.path(OUT, "package_versions.csv"))

tbl_path <- function(...) file.path(OUT, "tables", ...)
fig_path <- function(...) file.path(OUT, "figures", ...)
obj_path <- function(...) file.path(OUT, "objects", ...)

# Collapse an Ensembl-rowname matrix to HGNC symbols, summing duplicate symbols.
# Reads the stage-level `m` (objects/model.rds); every caller loads `m` at top
# level before calling this, so the lookup resolves in the global environment.
to_symbol <- function(mat) {
  s <- unname(m$symbols[rownames(mat)])
  keep <- !is.na(s) & s != "" & s != rownames(mat)
  rowsum(mat[keep, , drop = FALSE], group = s[keep])
}

# MSigDB access for Stage 03. It attaches msigdbr itself;
# the tryCatch covers the collection=/category= argument rename.
msig <- function(collection, sub = NULL) {
  x <- tryCatch(msigdbr(species = cfg$project$species, collection = collection,
                        subcollection = sub),
                error = function(e) msigdbr(species = cfg$project$species,
                                            category = collection, subcategory = sub))
  split(x$gene_symbol, x$gs_name)
}
flatten_le <- function(x) mutate(x, leadingEdge = vapply(leadingEdge, paste,
                                                        character(1), collapse = ";"))

# Build the five-state gene-expression profiles used consistently by Stages 07
# and 10. Input columns are individual reference cells and labels identify their
# GSE136831 lineage state. CPM is calculated per cell before state summaries.
make_GSE136831_profiles <- function(counts, labels) {
  counts <- as.matrix(counts)
  labels <- trimws(as.character(labels))
  expected <- unique(as.character(unlist(
    cfg$deconvolution$GSE136831_atlas$lineage_restricted_cell_types)))
  if (ncol(counts) != length(labels) || is.null(rownames(counts)))
    stop("GSE136831 reference counts/labels are malformed.", call. = FALSE)
  if (!setequal(unique(labels), expected))
    stop("GSE136831 reference states do not match config. ",
         set_diff_text(expected, unique(labels)), call. = FALSE)
  if (any(!is.finite(counts)) || any(counts < 0))
    stop("GSE136831 reference contains non-finite or negative counts.",
         call. = FALSE)
  lib <- colSums(counts)
  if (any(lib <= 0))
    stop("GSE136831 reference contains an empty cell.", call. = FALSE)
  cpm_cell <- sweep(counts, 2, lib, "/") * 1e6
  mean_cpm <- sapply(expected, function(state)
    rowMeans(cpm_cell[, labels == state, drop = FALSE]))
  detection <- sapply(expected, function(state)
    rowMeans(counts[, labels == state, drop = FALSE] > 0))
  dimnames(mean_cpm) <- list(rownames(counts), expected)
  dimnames(detection) <- list(rownames(counts), expected)
  list(
    mean_cpm = mean_cpm,
    detection_fraction = detection,
    state_cell_counts = setNames(
      as.integer(table(factor(labels, levels = expected))), expected),
    states = expected
  )
}

write_GSE136831_profiles <- function(profiles) {
  min_cpm <- cfg$crosstalk$GSE136831_min_cpm %||% 1
  min_detection <- cfg$crosstalk$GSE136831_min_detection_fraction %||% 0.10
  membership <- as.data.frame(profiles$mean_cpm) |>
    tibble::rownames_to_column("gene") |>
    tidyr::pivot_longer(-gene, names_to = "cell_type",
                        values_to = "mean_cpm") |>
    dplyr::mutate(
      detection_fraction = profiles$detection_fraction[
        cbind(gene, cell_type)],
      member = mean_cpm >= min_cpm &
        detection_fraction >= min_detection)
  profiles$membership <- membership
  profiles$thresholds <- list(
    min_cpm = min_cpm, min_detection_fraction = min_detection)
  reference_panel <- membership |>
    dplyr::filter(member) |>
    dplyr::mutate(
      n_reference_cells = unname(
        profiles$state_cell_counts[match(
          cell_type, names(profiles$state_cell_counts))])) |>
    dplyr::select(gene, cell_type, mean_cpm, detection_fraction,
                  n_reference_cells) |>
    dplyr::arrange(cell_type, dplyr::desc(mean_cpm), gene)
  saveRDS(profiles, obj_path("GSE136831_state_profiles.rds"))
  readr::write_csv(membership, tbl_path("GSE136831_gene_state_profiles.csv"))
  readr::write_csv(reference_panel,
                   tbl_path("GSE136831_reference_gene_panel.csv"))
  readr::write_csv(
    tibble::tibble(
      cell_type = names(profiles$state_cell_counts),
      n_cells = unname(profiles$state_cell_counts)),
    tbl_path("GSE136831_state_cell_counts.csv"))
  profiles
}

load_GSE136831_profiles <- function(rebuild = FALSE) {
  saved <- obj_path("GSE136831_state_profiles.rds")
  if (!rebuild && file.exists(saved)) {
    profiles <- readRDS(saved)
    wanted <- list(
      min_cpm = cfg$crosstalk$GSE136831_min_cpm %||% 1,
      min_detection_fraction =
        cfg$crosstalk$GSE136831_min_detection_fraction %||% 0.10)
    artifacts <- c(
      tbl_path("GSE136831_gene_state_profiles.csv"),
      tbl_path("GSE136831_reference_gene_panel.csv"),
      tbl_path("GSE136831_state_cell_counts.csv"))
    if (!identical(profiles$thresholds, wanted) ||
        any(!file.exists(artifacts)))
      profiles <- write_GSE136831_profiles(profiles)
    return(profiles)
  }
  configured <- cfg$crosstalk$GSE136831_reference_file %||% ""
  cap <- as.integer(cfg$deconvolution$GSE136831_atlas$source$
                      max_cells_per_donor_condition_cell_type)
  ref <- if (nzchar(configured)) configured else
    file.path(OUT, "cibersortx_upload",
              paste0("GSE136831_reference_donor_cap", cap, ".tsv"))
  if (!file.exists(ref))
    stop("GSE136831 single-cell reference is unavailable: ", ref,
         ". Run Stage 05 in this run or set ",
         "crosstalk$GSE136831_reference_file.", call. = FALSE)
  header <- strsplit(readLines(ref, n = 1), "\t", fixed = TRUE)[[1]]
  x <- read.delim(ref, skip = 1, header = FALSE, check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
  if (ncol(x) != length(header))
    stop("GSE136831 reference header/data column counts differ.", call. = FALSE)
  genes <- trimws(as.character(x[[1]]))
  if (any(!nzchar(genes)) || anyDuplicated(genes))
    stop("GSE136831 reference gene symbols are blank or duplicated.",
         call. = FALSE)
  counts <- as.matrix(x[-1])
  storage.mode(counts) <- "double"
  rownames(counts) <- genes
  write_GSE136831_profiles(make_GSE136831_profiles(counts, header[-1]))
}

# The config block is trusted here: enabled, dpi, and unique bare .svg
# filenames are set once in config.yml.
manuscript_png_config <- function() {
  png_cfg <- cfg$output$manuscript_png %||% list()
  list(
    enabled = isTRUE(png_cfg$enabled %||% FALSE),
    dpi = png_cfg$dpi %||% 600,
    subdir = png_cfg$subdir %||% "biorender_png",
    files = as.character(unlist(png_cfg$files %||% character(),
                                use.names = FALSE)))
}

publication_theme <- function() {
  theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    panel.background = element_rect(fill = "transparent", colour = NA),
    # Keep a crisp border around each plotting panel/facet while leaving the
    # outer canvas borderless and transparent for BioRender.
    panel.border = element_rect(
      fill = NA, colour = "grey55", linewidth = 0.3),
    # A separate axis line would double the left and bottom panel edges.
    axis.line = element_blank(),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.box.background = element_rect(fill = "transparent", colour = NA),
    strip.background = element_rect(
      fill = "grey95", colour = "grey55", linewidth = 0.3))
}

prepare_plot_for_export <- function(p) {
  if (inherits(p, "ggplot")) p + publication_theme() else p
}

with_figure_device <- function(open_device, draw) {
  open_device()
  tryCatch(draw(), finally = grDevices::dev.off())
  invisible(NULL)
}

# Save one transparent SVG and, for configured manuscript panels, a direct
# transparent PNG. `draw` is a zero-argument callback so grid/base graphics
# (pheatmap, circlize and composite figures) use the same export path.
save_figure <- function(file, width, height, draw) {
  if (!grepl("\\.svg$", file, ignore.case = TRUE))
    stop("Figure source filenames must end in .svg: ", file, call. = FALSE)
  svg_path <- fig_path(file)
  dir.create(dirname(svg_path), recursive = TRUE, showWarnings = FALSE)
  with_figure_device(
    function() svglite::svglite(svg_path, width = width, height = height,
                                bg = "transparent"),
    draw)

  png_cfg <- manuscript_png_config()
  if (png_cfg$enabled && basename(file) %in% png_cfg$files) {
    png_file <- paste0(tools::file_path_sans_ext(basename(file)), ".png")
    png_path <- fig_path(png_cfg$subdir, png_file)
    dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
    with_figure_device(
      function() grDevices::png(
        filename = png_path, width = width, height = height, units = "in",
        res = png_cfg$dpi, bg = "transparent"),
      draw)
  }
  invisible(svg_path)
}

save_plot <- function(p, file, width = 7, height = 5) {
  p <- prepare_plot_for_export(p)
  draw <- if (inherits(p, "pheatmap")) {
    function() {
      grid::grid.newpage()
      grid::grid.draw(p$gtable)
    }
  } else {
    function() print(p)
  }
  save_figure(file, width, height, draw)
  invisible(file)
}

# Colourblind-safe categorical palette (NPG, hues chosen to stay separable).
PAL <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")

# Symmetric blue-white-red fill for the Stage 03 Ryu leading-edge heatmap. The
# fallback limit keeps an all-zero or all-NA matrix from collapsing the scale.
diverging_fill <- function(x) {
  limit <- max(abs(x), na.rm = TRUE)
  if (!is.finite(limit) || limit == 0) limit <- 1
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                       midpoint = 0, limits = c(-limit, limit),
                       oob = scales::squish)
}

# Publication labels are deliberately separate from metadata/design IDs, so
# figure typography never changes a design-matrix column name or an
# external-file match.
condition_label_map <- function() {
  x <- unlist(cfg$display_labels$conditions %||% list(), use.names = TRUE)
  if (!length(x) || is.null(names(x)) || any(!nzchar(names(x))) ||
      any(!nzchar(as.character(x))))
    stop("display_labels$conditions must be a named, non-empty mapping.",
         call. = FALSE)
  setNames(as.character(x), names(x))
}

display_condition <- function(x) {
  map <- condition_label_map()
  raw <- as.character(x)
  labelled <- unname(map[raw])
  labelled[is.na(labelled)] <- raw[is.na(labelled)]
  labelled
}

display_condition_factor <- function(x) {
  map <- condition_label_map()
  factor(display_condition(x), levels = unname(map))
}

condition_palette <- function() {
  labels <- unname(condition_label_map())
  if (length(labels) > length(PAL))
    stop("Not enough colours for configured condition labels.", call. = FALSE)
  setNames(PAL[seq_along(labels)], labels)
}

comparison_label_map <- function() {
  values <- vapply(cfg$comparisons, function(cmp) {
    explicit <- cmp$label %||% ""
    if (nzchar(explicit)) return(as.character(explicit))
    paste(display_condition(cmp$numerator), "vs",
          display_condition(cmp$denominator))
  }, character(1))
  setNames(values,
           vapply(cfg$comparisons, function(cmp) cmp$name, character(1)))
}

display_comparison <- function(x) {
  map <- comparison_label_map()
  raw <- as.character(x)
  labelled <- unname(map[raw])
  labelled[is.na(labelled)] <- gsub("_", " ", raw[is.na(labelled)])
  labelled
}

reference_score_colours <- function(x) {
  if (!length(x) || any(!is.finite(x)))
    stop("Reference scores must be finite.", call. = FALSE)
  pal <- grDevices::colorRampPalette(c("#D9EAF7", "#67A9CF", "#2166AC"))(100)
  lim <- range(x)
  scaled <- if (diff(lim) == 0) rep(0.5, length(x)) else
    (x - lim[1]) / diff(lim)
  breaks <- seq(lim[1], lim[2], length.out = 3)
  list(
    colours = pal[pmax(1L, pmin(100L, 1L + floor(99 * scaled)))],
    legend = list(
      title = "reference score",
      labels = format(signif(breaks, 3), trim = TRUE),
      colours = pal[c(1, 50, 100)]))
}

#' Draw a ligand-receptor chord diagram with circlize.
#'
#' @param df data.frame(from, to, value); value sets the ribbon/segment width.
#' @param gene_ct named vector gene -> cell type (for colour and ordering).
#' @param ct_col named vector cell type -> colour.
#' @param file  file name inside the run's figures/ folder.
#' @param title figure title.
save_lr_flow <- function(df, gene_ct, ct_col, file, title = NULL,
                         width = 7.5, height = 7.5,
                         show_cell_type_prefix = TRUE,
                         inner_score = NULL) {
  if (!nrow(df)) return(invisible(NULL))
  if (!requireNamespace("circlize", quietly = TRUE))
    stop("Package `circlize` is required for ligand-receptor chord diagrams. ",
         "Run Rscript setup_packages.R first.", call. = FALSE)
  if (!is.null(inner_score) &&
      (length(inner_score) != nrow(df) ||
       any(!is.finite(inner_score)) || any(inner_score < 0)))
    stop("inner_score must contain one finite, non-negative reference score ",
         "per L-R row.", call. = FALSE)
  ct_levels <- names(ct_col)
  display_label <- function(x) {
    if (show_cell_type_prefix) x else sub("^[^:]+:", "", x)
  }
  # Preserve the original layout: group both senders and receivers by state.
  sec_order <- names(gene_ct)[order(match(gene_ct, ct_levels), names(gene_ct))]
  grid_col <- setNames(unname(ct_col[gene_ct]), names(gene_ct))
  sender_ct <- unname(gene_ct[df$from])
  if (anyNA(sender_ct) || any(!sender_ct %in% names(ct_col)))
    stop("Every chord sender must have a configured cell-type colour.",
         call. = FALSE)
  sender_col <- grDevices::adjustcolor(
    unname(ct_col[sender_ct]), alpha.f = 0.62)
  score_style <- if (is.null(inner_score)) NULL else
    reference_score_colours(inner_score)
  cex <- max(0.4, min(0.8, 24 / length(sec_order)))
  draw <- function() {
    circlize::circos.clear()
    on.exit(circlize::circos.clear(), add = TRUE)
    circlize::circos.par(points.overflow.warning = FALSE)
    chord_pos <- circlize::chordDiagram(
      df, grid.col = grid_col, order = sec_order, directional = 1,
      direction.type = c("diffHeight", "arrows"),
      link.arr.type = "triangle", link.arr.length = 0.12,
      link.arr.width = 0.06, link.arr.col = sender_col,
      link.arr.lwd = 0, link.border = NA,
      annotationTrack = "grid", link.sort = FALSE,
      preAllocateTracks = list(track.height = 0.12),
      col = sender_col)
    if (!is.null(score_style)) {
      if (nrow(chord_pos) != nrow(df) ||
          !identical(as.character(chord_pos$rn), as.character(df$from)) ||
          !identical(as.character(chord_pos$cn), as.character(df$to)))
        stop("circlize changed L-R row order; refusing to misalign the ",
             "pair-level reference-score ring.", call. = FALSE)
      inner_track <- NA_integer_
      circlize::circos.track(
        ylim = c(0, 1), track.height = 0.035,
        cell.padding = c(0, 0, 0, 0), bg.border = NA,
        panel.fun = function(x, y)
          inner_track <<- circlize::CELL_META$track.index)
      for (i in seq_len(nrow(chord_pos))) {
        endpoint_col <- score_style$colours[i]
        circlize::circos.rect(
          chord_pos$x1[i] - abs(chord_pos$value1[i]), 0,
          chord_pos$x1[i], 1,
          col = endpoint_col, border = endpoint_col,
          sector.index = chord_pos$rn[i], track.index = inner_track)
        circlize::circos.rect(
          chord_pos$x2[i] - abs(chord_pos$value2[i]), 0,
          chord_pos$x2[i], 1,
          col = endpoint_col, border = endpoint_col,
          sector.index = chord_pos$cn[i], track.index = inner_track)
      }
    }
    circlize::circos.track(track.index = 1, bg.border = NA,
                           panel.fun = function(x, y)
      circlize::circos.text(circlize::CELL_META$xcenter,
                            circlize::CELL_META$ylim[1],
                            display_label(circlize::CELL_META$sector.index),
                            facing = "clockwise", niceFacing = TRUE,
                            adj = c(0, 0.5), cex = cex))
    graphics::legend("bottomleft", legend = names(ct_col), fill = ct_col,
                     bty = "n", cex = 0.8, title = "cell type")
    if (!is.null(score_style))
      graphics::legend(
        "bottomright", legend = score_style$legend$labels,
        fill = score_style$legend$colours, bty = "n", cex = 0.75,
        title = score_style$legend$title)
    if (!is.null(title)) graphics::title(title, cex.main = 0.95)
  }
  save_figure(file, width, height, draw)
  invisible(fig_path(file))
}
