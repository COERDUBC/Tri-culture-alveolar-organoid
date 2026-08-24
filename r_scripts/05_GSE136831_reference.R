#!/usr/bin/env Rscript
# 05 — Build the GSE136831 single-cell reference for CIBERSORTx from the four
# processed GSE136831 files (Adams et al. 2020). Run once, and again whenever
# `lineage_restricted_cell_types` changes.
#
# Restricted to the states the tri-culture can contain plus the IPF states the
# model should reproduce; airway states are excluded by experimental design, so
# they cannot compete for signal in the fit. Sampling is capped per donor per
# cell type so no single donor dominates a state; one reference is written per
# configured cap, the first being the primary and the rest the reference-size
# sensitivity analysis.
#
# Reads the 2 GB compressed matrix: run in a session with >= 64 GB RAM.

source("r_scripts/common.R")
suppressPackageStartupMessages({ library(Matrix); library(tidyr) })

kcfg <- cfg$deconvolution$GSE136831_atlas
scfg <- kcfg$source
up <- file.path(OUT, "cibersortx_upload")
dir.create(up, recursive = TRUE, showWarnings = FALSE)

paths <- c(matrix = scfg$matrix_file, genes = scfg$gene_file,
           barcodes = scfg$barcode_file, metadata = scfg$metadata_file)
paths <- setNames(file.path(scfg$data_dir, unname(paths)), names(paths))
if (any(!file.exists(paths)))
  stop("Missing GSE136831 file(s): ", paste(paths[!file.exists(paths)], collapse = ", "))

gene_info <- read_tsv(paths[["genes"]], show_col_types = FALSE, progress = FALSE)
meta <- read_tsv(paths[["metadata"]], show_col_types = FALSE, progress = FALSE)
barcodes <- read_lines(paths[["barcodes"]], progress = FALSE)

gene_col <- scfg$gene_symbol_column
cell_id_col <- scfg$barcode_column
needed <- c(cell_id_col, scfg$cell_type_column, scfg$donor_column,
            scfg$condition_column)
if (!gene_col %in% names(gene_info))
  stop("Gene file lacks ", gene_col, ". Available: ",
       paste(names(gene_info), collapse = ", "))
if (any(!needed %in% names(meta)))
  stop("Metadata lacks: ", paste(setdiff(needed, names(meta)), collapse = ", "))

# Some exports carry a one-line barcode header.
if (length(barcodes) == nrow(meta) + 1L &&
    tolower(gsub("[^[:alnum:]]", "", barcodes[1])) %in%
      c("cellbarcode", "cellbarcodeidentity")) barcodes <- barcodes[-1]
ord <- match(barcodes, meta[[cell_id_col]])
if (length(barcodes) != nrow(meta) || anyNA(ord) || anyDuplicated(barcodes))
  stop("Barcode and metadata files do not describe the same cells.")

meta <- meta[ord, , drop = FALSE]
meta$matrix_col <- seq_along(barcodes)
meta$cell_id <- barcodes
meta$cell_type_raw <- trimws(as.character(meta[[scfg$cell_type_column]]))
meta$cell_type <- gsub("_+", "_", gsub("[^[:alnum:]_]+", "_", meta$cell_type_raw))
meta$donor <- trimws(as.character(meta[[scfg$donor_column]]))
meta$condition <- trimws(as.character(meta[[scfg$condition_column]]))
meta <- meta |>
  filter(!is.na(cell_type), cell_type != "", !is.na(donor), donor != "",
         condition %in% unlist(scfg$include_conditions),
         !cell_type_raw %in% trimws(as.character(unlist(scfg$exclude_cell_types %||% character()))))

# A wrong Disease_Identity literal would silently empty `meta` and then surface
# as a confusing "no cell types available" error further down.
missing_conditions <- setdiff(
  trimws(as.character(unlist(scfg$include_conditions))), unique(meta$condition))
if (length(missing_conditions))
  stop("Configured condition(s) absent from the atlas: ",
       paste(missing_conditions, collapse = ", "), "\nAvailable: ",
       paste(sort(unique(meta$condition)), collapse = ", "))

# All eligible donors, before restricting to the five selected states
all_donors <- sort(unique(meta$donor))
donor_condition <- distinct(meta, donor, condition)

# Keep only the configured states, and fail loudly if one is absent.
wanted <- unique(trimws(as.character(unlist(kcfg$lineage_restricted_cell_types))))
missing <- setdiff(wanted, unique(meta$cell_type))
if (length(missing))
  stop("Configured cell type(s) absent from the atlas: ",
       paste(missing, collapse = ", "), "\nAvailable: ",
       paste(sort(unique(meta$cell_type)), collapse = ", "))
meta <- filter(meta, cell_type %in% wanted)

# -------------------------------------------------------------------------
# Inspect donor representation before downsampling
# -------------------------------------------------------------------------

donor_state_counts <- meta |>
  count(donor, cell_type, name = "n_cells")

# Cell type × donor cell-count table
donor_state_matrix <- tidyr::expand_grid(
  donor = all_donors,
  cell_type = wanted
) |>
  left_join(
    donor_state_counts,
    by = c("donor", "cell_type")
  ) |>
  mutate(n_cells = replace_na(n_cells, 0L)) |>
  pivot_wider(
    names_from = cell_type,
    values_from = n_cells
  ) |>
  left_join(donor_condition, by = "donor") |>
  relocate(condition, .after = donor) |>
  arrange(condition, donor)

write_csv(
  donor_state_matrix,
  file.path(up, "celltype_by_donor_counts.csv")
)
message(nrow(donor_state_matrix), " donors x ", length(wanted),
        " states -> celltype_by_donor_counts.csv")
# full table with 0
donor_state_complete <- tidyr::expand_grid(
  donor = all_donors,
  cell_type = wanted
) |>
  left_join(
    donor_state_counts,
    by = c("donor", "cell_type")
  ) |>
  mutate(n_cells = replace_na(n_cells, 0L))

# summary table for the scRNA reference data
donor_state_summary <- donor_state_complete |>
  group_by(cell_type) |>
  summarise(
    total_donors = n_distinct(donor),
    donors_with_state = sum(n_cells > 0),
    donors_without_state = sum(n_cells == 0),
    donor_coverage_pct = 100 * mean(n_cells > 0),
    total_cells = sum(n_cells),
    min_cells_per_donor = min(n_cells[n_cells > 0]),
    median_cells_per_donor = median(n_cells[n_cells > 0]),
    mean_cells_per_donor = mean(n_cells[n_cells > 0]),
    max_cells_per_donor = max(n_cells),
    .groups = "drop"
  )

write_csv(
  donor_state_summary,
  file.path(up, "donor_state_summary.csv")
)
message("Per-state donor coverage -> donor_state_summary.csv")

min_cells <- as.integer(scfg$min_cells_per_cell_type %||% 3)
available <- count(meta, cell_type, condition, name = "available_cells")
small <- count(meta, cell_type, name = "n") |> filter(n < min_cells)
if (nrow(small))
  stop("Cell type(s) with < ", min_cells, " cells: ",
       paste(small$cell_type, collapse = ", "))

# --- sampling ---------------------------------------------------------------
# Draw up to `cap` cells per donor per cell type, so every donor carries equal
# weight and the reference describes the donor population instead of the few
# donors that happen to contribute the most cells of a state. `condition` is not
# a grouping variable — a donor belongs to exactly one condition, so adding it
# would not change the groups.
cap <- as.integer(scfg$max_cells_per_donor_condition_cell_type)
set.seed(SEED)
selected <- meta |>
  group_by(donor, cell_type) |>
  slice_sample(n = cap) |>
  ungroup() |>
  arrange(cell_type, condition, donor, cell_id)
if (length(setdiff(wanted, unique(selected$cell_type))))
  stop("Balanced sampling dropped a configured cell type.")
message("donor cap ", cap, ": ", nrow(selected), " cells across ",
        n_distinct(selected$cell_type), " cell types")
print(count(selected, cell_type, condition, name = "n_cells"))

# --- counts ----------------------------------------------------------------
message("Reading the compressed GSE136831 matrix (needs >= 64 GB RAM)...")
con <- gzfile(paths[["matrix"]], "rt")
counts <- tryCatch(readMM(con), finally = close(con))
if (nrow(counts) != nrow(gene_info) || ncol(counts) != length(barcodes))
  stop("Matrix dims ", nrow(counts), "x", ncol(counts), " do not match genes (",
       nrow(gene_info), ") / barcodes (", length(barcodes), ").")
counts <- as(counts[, selected$matrix_col, drop = FALSE], "CsparseMatrix")
invisible(gc())

# CIBERSORTx keys on gene symbols: drop unnamed genes, sum duplicated symbols.
sym <- trimws(as.character(gene_info[[gene_col]]))
keep <- !is.na(sym) & sym != ""
counts <- counts[keep, , drop = FALSE]; sym <- sym[keep]
uniq <- unique(sym)
counts <- sparseMatrix(i = match(sym, uniq), j = seq_along(sym), x = 1,
                       dims = c(length(uniq), length(sym))) %*% counts
rownames(counts) <- uniq
colnames(counts) <- selected$cell_id
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

# --- write the reference ---------------------------------------------------
# Header is the cell-type label per column (what CIBERSORTx groups on). Written
# in row blocks so the dense conversion never holds the whole matrix.
suffix <- paste0("_donor_cap", cap)
ref_path <- file.path(up, paste0("GSE136831_reference", suffix, ".tsv"))
con <- file(ref_path, "wt")
tryCatch({
  writeLines(paste(c("GeneSymbol", selected$cell_type), collapse = "\t"), con)
  rows <- seq_len(nrow(counts))
  for (first in seq(1L, length(rows), by = 250L)) {
    blk <- rows[first:min(first + 249L, length(rows))]
    write.table(as.matrix(counts[blk, , drop = FALSE]), con, sep = "\t",
                quote = FALSE, row.names = TRUE, col.names = FALSE)
  }
}, finally = close(con))

write_csv(select(selected, cell_id, cell_type, cell_type_raw, donor, condition),
          file.path(up, paste0("GSE136831_manifest", suffix, ".csv")))
write_csv(full_join(available, count(selected, cell_type, condition,
                                     name = "selected_cells"),
                    by = c("cell_type", "condition")) |>
            mutate(selected_cells = coalesce(selected_cells, 0L)) |>
            arrange(cell_type, condition),
          file.path(up, paste0("GSE136831_sampling_summary", suffix, ".csv")))

# State-level expression and detection profiles for the downstream
# connectomeDB/NicheNet analyses. invisible() because this returns the whole
# profiles list, and a top-level visible value would be auto-printed by Rscript.
invisible(write_GSE136831_profiles(
  make_GSE136831_profiles(as.matrix(counts), selected$cell_type)))

cat(sprintf('
Wrote %s
  %d genes x %d cells, %d cell types, donor cap %d

Upload to https://cibersortx.stanford.edu as a "Single Cell Reference Matrix",
run Impute Cell Fractions against bulk_mixture_cpm.tsv (from Stage 04) with
permutations = %s, quantile normalisation = %s, batch correction = %s,
then set deconvolution$cibersortx$GSE136831_lineage_result in config.yml.
',
  ref_path, nrow(counts), ncol(counts), n_distinct(selected$cell_type), cap,
  cfg$deconvolution$permutations %||% 1000,
  if (isTRUE(cfg$deconvolution$quantile_normalization)) "on" else "off",
  cfg$deconvolution$cibersortx$GSE136831_batch_correction %||% "S-mode"))

message("05 done -> ", up)
