## data-raw/build_htr_scores.R
## -----------------------------------------------------------------------------
## Build a per-cell HTR receptor expression table from the macaque PatchSeq
## Seurat object. Output is the data foundation for manuscript figures 2, 3, 4,
## 6, and S2.
##
## Run from the proj5HT repo root:
##   Rscript data-raw/build_htr_scores.R
##
## Output:
##   data-raw/htr_scores.csv         per-cell expression + family sums
##   data-raw/htr_scores.rds         same, as a data.frame
##
## Scaffolding notes (2026-06-01):
##   - HTR1A and HTR1B values in this build come straight from the primary
##     macaque PatchSeq Seurat object and are KNOWN-INACCURATE; Scott will
##     overlay corrected values from accessory files in a follow-up pass.
##   - Two columns flag this: `htr1a_status` and `htr1b_status` are set to
##     "provisional_seurat_primary". After patching, set them to the
##     accessory-file tag (e.g. "patched_<source>_<date>").
##   - All other HTR genes (HTR1D/E/F, HTR2A/B/C, HTR4/5A/6/7) are treated
##     as final from the primary Seurat object.
##   - Subclass identity comes from the headBCI hierarchical classifier
##     (`hier_subclass_assignment` in `headBCI::sheets$cleanMD`), attached
##     via `proj5HT::add_hier_to_seurat()`. The Seurat object's legacy
##     `subclass_Corr` column is preserved but no longer canonical.
##   - Expression units: the Seurat object's `counts` layer is curated as
##     counts-per-million (CPM), NOT raw counts. Do NOT run
##     Seurat::NormalizeData() on it — values are already on a comparable
##     per-cell scale. FetchData() emits a harmless 'data layer not found,
##     counts used' warning which is the desired behaviour here.
##   - Family sums and HTR_balance are recomputed every run from the current
##     values, so re-running this script after Scott patches HTR1A/HTR1B
##     will refresh the derived columns automatically.
## -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(headBCI)
  ## load proj5HT in-place so this script picks up uncommitted source
  ## changes (e.g. to R/hier_metadata5HT.R) without needing devtools::install().
  devtools::load_all(".", quiet = TRUE, export_all = FALSE)
})

## ---- config ----------------------------------------------------------------

SEURAT_PATH <- "data-raw/Seurat/offPipeline-Macaque_PatchSeq.rds"
OUT_CSV     <- "data-raw/htr_scores.csv"
OUT_RDS     <- "data-raw/htr_scores.rds"

## Receptor groupings — must match manu/_shared/setup.R HTR_FAMILIES.
HTR1 <- c("HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F")
HTR2 <- c("HTR2A", "HTR2B", "HTR2C")
HTR0 <- c("HTR4",  "HTR5A", "HTR6",  "HTR7")
HTR_ALL <- c(HTR1, HTR2, HTR0)

## Genes whose values we KNOW are inaccurate in the primary Seurat object
## and will be patched from accessory files later. Tagged in `*_status` cols.
PROVISIONAL_GENES <- c("HTR1A", "HTR1B")
PROVISIONAL_TAG   <- "provisional_seurat_primary"

## ---- load Seurat object ----------------------------------------------------

message("Loading Seurat object: ", SEURAT_PATH)
obj <- readRDS(SEURAT_PATH)

## attach hier_* columns from headBCI::sheets$cleanMD (canonical subclass)
obj <- proj5HT::add_hier_to_seurat(obj)

## sanity check — which HTR genes are present in the assay?
present <- intersect(HTR_ALL, rownames(obj))
missing <- setdiff(HTR_ALL, rownames(obj))
if (length(missing)) {
  warning("HTR genes missing from Seurat object (filled with NA): ",
          paste(missing, collapse = ", "))
}

## ---- per-cell expression ---------------------------------------------------

## Pull normalized expression for all HTR genes plus the metadata columns
## the manuscript scripts use as cell-level keys.
meta_cols <- c("cell_name",
               ## canonical subclass identity (headBCI hier classifier)
               "hier_class_assignment",
               "hier_neighborhood_assignment",
               "hier_subclass_assignment",
               "hier_cluster_assignment",
               "hier_subclass_bootstrap_prob",
               "hier_subclass_avg_corr",
               ## legacy Seurat-internal calls (kept for QC comparison)
               "subclass_Corr", "subclass_Tree", "cluster_Corr_label",
               "predicted_subclass",
               ## depth + cortical area
               "Cortical_area", "Species",
               "LIMS_depth", "Depth_from_pia",
               ## QC
               "Genes.Detected", "percent_reads_aligned_to_introns",
               "marker_sum_norm_label", "quality_score_label",
               "score.Corr_label")
meta_cols <- intersect(meta_cols, colnames(obj@meta.data))

df <- Seurat::FetchData(obj, vars = c(meta_cols, present))

## normalize hier assignment labels (single source of truth)
df <- proj5HT::normalize_hier_assignments(df)

## fill any genes missing from the object with NA so downstream code is stable
for (g in missing) df[[g]] <- NA_real_

## reorder columns: metadata first, then HTR genes in canonical order
df <- df[, c(meta_cols, HTR_ALL)]

## ---- provisional-source flags ----------------------------------------------

for (g in PROVISIONAL_GENES) {
  status_col <- paste0(tolower(g), "_status")
  df[[status_col]] <- PROVISIONAL_TAG
}

## ---- derived family sums + balance score -----------------------------------

## NA-tolerant rowSums: treat NA as 0 only if at least one member is observed;
## otherwise leave NA so downstream callers can detect missing data.
family_sum <- function(df, genes) {
  m <- as.matrix(df[, genes, drop = FALSE])
  has_any <- rowSums(!is.na(m)) > 0
  out <- rowSums(m, na.rm = TRUE)
  out[!has_any] <- NA_real_
  out
}

df$HTR1_sum    <- family_sum(df, HTR1)
df$HTR2_sum    <- family_sum(df, HTR2)
df$HTR0_sum    <- family_sum(df, HTR0)
df$HTR_balance <- df$HTR2_sum - df$HTR1_sum   ## > 0 → excitatory bias

## ---- write outputs ---------------------------------------------------------

if (!dir.exists("data-raw")) dir.create("data-raw")

readr::write_csv(df, OUT_CSV)
saveRDS(df, OUT_RDS)

message(sprintf(
  "Wrote %d cells x %d cols\n  csv: %s\n  rds: %s",
  nrow(df), ncol(df), OUT_CSV, OUT_RDS
))
message("Provisional columns (need accessory-file patch): ",
        paste(PROVISIONAL_GENES, collapse = ", "))
