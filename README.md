# Tri-culture alveolar-fibrosis organoid — RNA-seq analysis

Analysis code for the manuscript of Characterization of tri-culture alveolar-like organoids to study fibrotic lung diseases. Read-level processing and quantification
(nf-core/rnaseq 3.26.0, STAR 2.7.11b, Salmon 1.10.3, GRCh38 / GENCODE 49) run
upstream; every stage here starts from the completed STAR-Salmon output.

A tri-culture alveolar organoid (one fibroblast and two epithelial cell lines)
treated with **TGF-β1** or **Bleomycin** against untreated **Control**.

## Design

18 libraries: 3 conditions × 3 replicates × 2 batches. `metadata.csv` carries
`sample, condition, batch` with conditions `Control`, `TGF`, `Bleo`.

The design is `~ batch + condition`, so batch is a covariate in every test;
`removeBatchEffect` is applied only to the log-CPM used for the PCA and heatmaps. 
`Control` is the reference level. Inference is limma-voom with eBayes moderation, 
and GSEA ranks genes by the same moderated t-statistic. 
Contrasts: TGF-β1 vs Control, Bleomycin vs Control, TGF-β1 vs Bleomycin.

Figures display the internal IDs as TGF-β1 / Bleomycin / Control through
`display_labels` in `config.yml`; do not rename the metadata IDs.

## Layout

```
config.yml           # committed template: paths, design, thresholds, external references
config_local.yml     # same schema with this machine's real paths — used if present, git-ignored
setup_packages.R     # one-time package installation
r_scripts/
  common.R           # config, run folder, save_plot, palettes, GSE136831 state profiles
  01_import.R        # tximport + limma-voom fit; PCA, sample distances
  02_de.R            # per-contrast DE, volcano + heatmap, targeted marker panel
  03_gsea.R          # MSigDB GSEA + the Ryu 2024 aberrant-basaloid module
  04_cibersortx_prep.R      # bulk mixture file for the CIBERSORTx website
  05_GSE136831_reference.R  # donor-balanced five-state reference (>= 64 GB RAM)
  06_deconvolution.R # five-state abundance + treatment tests
  07_lr.R            # connectomeDB2025 ligand-receptor candidates
  08_nichenet.R      # NicheNet ligand activity + active-LR chord
bash_run/            # step1 / step2 / step3 — plain bash, no scheduler
data/                # CIBERSORTx signature matrix used for the reported deconvolution
session_info/        # R and package versions from the run that produced the results
runs/                # per-run pipeline output (git-ignored)
```

## Setup

``` bash
Rscript setup_packages.R
```

Then set the input and external-reference paths. `config.yml` is the committed
template and carries `<ProjectPath>` placeholders; copy it to `config_local.yml`
and put the real paths there — that file is git-ignored and every stage uses it
in preference to `config.yml`, so local paths never show up in a diff. Set
`CONFIG=<file>` to force a specific one. Each stage prints the file it used and
copies it into the run folder as `config.yml`. The pipeline
reads the Salmon `quant.sf` files, `salmon.merged.tx2gene.tsv`, and the metadata
CSV. No genome FASTA or GTF is used. Gene symbols come from the tx2gene third
column, then the offline OrgDb, then the Ensembl ID.

## Run

From the repository root (RStudio Terminal or any shell). Export the UTF-8
locale first — without it svglite renders `TGF-β1` as mojibake:

``` bash
cd /path/to/Tri_organoid
export LANG=C.UTF-8 LC_ALL=C.UTF-8

# step 1 — stages 01-03: DE and pathways
Rscript r_scripts/01_import.R            # creates runs/run_<timestamp>/, writes runs/CURRENT_RUN
Rscript r_scripts/02_de.R                # per-contrast DE, volcanoes, heatmaps, marker panel
Rscript r_scripts/03_gsea.R              # MSigDB GSEA + Ryu 2024 aberrant-basaloid module

# step 2 — stages 04-05: CIBERSORTx inputs, >= 64 GB RAM for 05
Rscript r_scripts/04_cibersortx_prep.R   # bulk_mixture_cpm.tsv
Rscript r_scripts/05_GSE136831_reference.R  # GSE136831_reference_donor_cap<k>.tsv

# ... run the two CIBERSORTx jobs, then set
# deconvolution.cibersortx.GSE136831_lineage_result in config_local.yml ...

# step 3 — stages 06-08: deconvolution and cross-talk
Rscript r_scripts/06_deconvolution.R     # five-state abundance + treatment tests
Rscript r_scripts/07_lr.R                # connectomeDB2025 ligand-receptor candidates
Rscript r_scripts/08_nichenet.R          # NicheNet ligand activity + active-LR chord
```

The three wrapper scripts are equivalent — each sets the locale itself and calls
its stages in order:

``` bash
bash bash_run/step1_de.sh                 # stages 01-03: DE and pathways
bash bash_run/step2_cibersortx_prep.sh    # stages 04-05: CIBERSORTx inputs, >= 64 GB RAM
# ... run the CIBERSORTx job, set its result path in config_local.yml ...
bash bash_run/step3_downstream.sh         # stages 06-08: deconvolution and cross-talk
```

`r_scripts/common.R` is sourced by every stage, never run on its own. Stage 01
creates the run folder and writes `runs/CURRENT_RUN`; every later stage reads
that pointer, so run 01 first, and keep the stages in numerical order within a step. 
Confirm stage 01 prints `Model: NNNNN genes x 18 samples`. 
Stages 02-08 can be re-run individually against the current run folder 
(for example after changing a threshold in `config.yml`) 
without repeating stage 01; to target an older run instead of the
one in `runs/CURRENT_RUN`, pass its folder name in `OUT`, e.g.
`OUT=run_2026-08-20_10-00-00 Rscript r_scripts/03_gsea.R`. Stage 08 prints a
skip message instead of failing when its NicheNet files are blank or still
`TO_FILL`.

Each run folder holds `objects/` (RDS cache), `tables/`, `figures/`
(+ `figures/biorender_png/`), the `config.yml` used, `package_versions.csv` and
`sessionInfo.txt`. To keep a console log: `bash bash_run/step1_de.sh 2>&1 | tee step1.log`.
The copies in `session_info/` are from the run that produced the reported
results.

### CIBERSORTx

CIBERSORTx is a hosted tool, so it runs between steps 2 and 3. Step 2 writes
`runs/<run>/cibersortx_upload/` containing the mixture `bulk_mixture_cpm.tsv`
(stage 04) and the single-cell reference `GSE136831_reference_donor_cap<k>.tsv`
(stage 05, where `<k>` is the configured cells-per-donor cap). At
https://cibersortx.stanford.edu run two jobs: **Create Signature Matrix**
(Custom / scRNA-seq, that reference, quantile normalization off), then **Impute
Cell Fractions** (that signature, the mixture, S-mode batch correction with the
same reference, 1000 permutations, quantile normalization off, absolute mode on).
Absolute mode is deliberate — relative fractions are compositional, so one state
can only rise if others fall. Then set
`deconvolution.cibersortx.GSE136831_lineage_result` (and
`crosstalk.GSE136831_reference_file`) in `config_local.yml` and run step 3. These
settings mirror `deconvolution:` in `config.yml`; change them in both places.
Step 2 can be skipped on a re-run once both files exist and are configured.

The signature matrix behind the reported results is committed at
`data/Triple_coculture_3D_Alveolar_Organoid_Deconvolution_SigMatrix.txt`
(2,210 genes × the five states), so the imputation can be repeated without
rebuilding the single-cell reference or re-running Create Signature Matrix.

## Figures

All cell-state analyses use five states: ATI, ATII, Aberrant_Basaloid,
Fibroblast, Myofibroblast.

| Figure | Script | Table |
|---|---|---|
| `pca.svg`, `sample_distances.svg` | `01_import.R` | — |
| `volcano_{TGF,Bleo}_vs_Control.svg`, `heatmap_{TGF,Bleo}_vs_Control.svg` | `02_de.R` | `de_<contrast>.csv`, `de_significant_all_contrasts.csv` |
| `targeted_aberrant_basal_marker_expression.svg` | `02_de.R` | `targeted_marker_expression{,_summary}.csv`, `targeted_marker_treatment_effects.csv` |
| `gsea_hallmark_dotplot_{TGF,Bleo}_vs_Control.svg` | `03_gsea.R` | `gsea_msigdb.csv` |
| `gsea_ryu2024_programme_and_leading_edge.svg` | `03_gsea.R` | `gsea_ryu2024.csv`, `ryu2024_leading_edge{,_heatmap}.csv` |
| `deconvolution_GSE136831_states{,_boxplot}.svg` | `06_deconvolution.R` | `deconvolution_GSE136831_states.csv`, `deconvolution_treatment_tests.csv` |
| `lr_treatment_associated.svg` | `07_lr.R` | `lr_candidates.csv` |
| `nichenet_lr_chord_{TGF,Bleo}_vs_Control.svg` | `08_nichenet.R` | `nichenet_ligand_activity_<contrast>.csv`, `nichenet_active_lr_<contrast>.csv` |

## External resources

Downloaded once and pinned in `config.yml`; nothing is fetched at run time.

| Resource | Used by | Source |
|---|---|---|
| GSE136831 IPF Cell Atlas (4 processed files) | `05`-`08` | GEO GSE136831 |
| connectomeDB2025 human ligand-receptor table | `07` | https://connectomedb.org |
| NicheNet prior model (ligand-target matrix, LR network) | `08` | nichenetr Zenodo record |
| Ryu et al. 2024 aberrant-basaloid programme (18-gene subset, listed in `config.yml`) | `03` | doi:10.1186/s12864-024-11031-5 |

## Data availability

The sequencing data supporting this analysis are available from the
corresponding author on reasonable request.
