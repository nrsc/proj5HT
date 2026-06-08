# ============================================================================
# dev/Seurat/Seurat_Objects.R
# ----------------------------------------------------------------------------
# Build Seurat objects for the proj5HT cohorts (Patch-Seq + 10x FACS,
# macaque / human / chimp / gorilla, on- and off-pipeline subsets).
#
# Refactor (2026-06-02):
#   - Collapses the 8 near-identical SeuratHCT_* functions in the
#     pre-refactor version (see
#     .archive/2026-06-02_ketwash_spikepuff/Seurat_Objects.pre-refactor.R)
#     into a single internal core, `build_seurat_from_feather()`.
#   - Per-cohort builders are thin wrappers (~10 lines each).
#   - `add_hier_to_seurat()` is invoked automatically when the
#     annotation feather carries `cell_name_label` (i.e. Patch-Seq
#     cohorts), so hier_* columns land in `@meta.data` at build time.
#   - All cohorts are saved under `data-raw/Seurat/` with a single tag
#     scheme: `<species>_<modality>[_<subset>].rds`.
#   - The auto-execution at source-time is removed; use
#     `build_all_seurat_objects()` (or call individual builders) once
#     you actually want to rebuild.
#
# Conventions:
#   - All counts feathers are assumed to be cell-id-columned tables with
#     a leading `gene` column; CPM normalization + log2(x+1) is applied.
#   - Patch-Seq cohorts join against `headBCI::sheets$cleanMD` (replaces
#     the previous `projHCT::sheets$MD` dependency).
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(feather)
  library(dplyr)
})

## ----------------------------------------------------------------------------
## Default cohort knobs
## ----------------------------------------------------------------------------

.SEU_PATHS <- list(
  macaque_patchseq = "~/proj5HT/den/macaque/MTG/GreatApes_Macaque_NCBI_RSC-204-400_map_full/",
  human_patchseq   = "~/proj5HT/den/macaque/MTG/GreatApes_Human_RSC-122-380_map_full/",
  macaque_10x      = "~/proj5HT/den/GreatApes_Macaque_NCBI/",
  human_10x        = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Human/",
  chimp_10x        = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Chimp/",
  gorilla_10x      = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Gorilla/"
)

.SEU_OUTDIR <- "~/proj5HT/data-raw/Seurat"

.EXC_SUBCLASS <- c("L5 IT", "L5 ET", "L5/6 NP", "L4 IT",
                   "L2/3 IT", "L6 IT", "L6 CT", "L6b")
.EXC_INH_SUBCLASS <- c(.EXC_SUBCLASS, "Chandelier", "Pvalb",
                       "Vip", "Sst", "Sst Chodl")
.EXC_NEIGHBORHOODS <- c("it_types", "l5et_l56np_l6ct_l6b")

## ----------------------------------------------------------------------------
## Core builder
## ----------------------------------------------------------------------------

#' Build one Seurat object from a `{counts, anno}` feather pair.
#'
#' Cell-name detection:
#'   - Patch-Seq feathers carry `cell_name_label` → triggers the join to
#'     `headBCI::sheets$cleanMD` and a follow-on `add_hier_to_seurat()`.
#'   - 10x FACS feathers don't carry per-cell ephys metadata; the join
#'     is skipped and only the feather's own annotation is attached.
#'
#' @param path Directory containing `counts.feather` (or `data.feather`)
#'   and `anno.feather`.
#' @param subclass Optional character vector of subclass labels to keep.
#'   Matches against `subclass_label` (10x) or `subclass_Corr_label`
#'   (Patch-Seq), whichever is present.
#' @param neighborhood Optional character vector of neighborhood labels
#'   (matched against `neighborhood_label` / `neighborhood_Corr_label`).
#' @param counts_file Filename of the counts feather (default tries
#'   `counts.feather` then `data.feather`).
#' @param normalize One of `"cpm_log2"` (10x raw counts → CPM → log2)
#'   or `"log2_only"` (Patch-Seq feathers that are already RPM).
#' @param project Seurat project name.
#' @return Seurat object. Patch-Seq cohorts (annotation feather carrying
#'   `cell_name_label`) automatically get `headBCI::sheets$cleanMD` joined
#'   in — hier_* columns included — with assignment labels normalized via
#'   `proj5HT::normalize_hier_assignments()`.
#' @keywords internal
build_seurat_from_feather <- function(path,
                                      subclass     = NULL,
                                      neighborhood = NULL,
                                      counts_file  = NULL,
                                      normalize    = c("cpm_log2", "log2_only"),
                                      project      = "proj5HT_seurat") {

  normalize <- match.arg(normalize)
  stopifnot(dir.exists(path))

  ## ---- annotation ----------------------------------------------------------
  anno <- feather::read_feather(file.path(path, "anno.feather"))
  anno <- as.data.frame(anno)

  sub_col <- if ("subclass_Corr_label" %in% names(anno)) "subclass_Corr_label" else "subclass_label"
  nbh_col <- if ("neighborhood_Corr_label" %in% names(anno)) "neighborhood_Corr_label" else "neighborhood_label"

  if (!is.null(subclass) && sub_col %in% names(anno)) {
    anno <- anno[anno[[sub_col]] %in% subclass, , drop = FALSE]
  }
  if (!is.null(neighborhood) && nbh_col %in% names(anno)) {
    anno <- anno[anno[[nbh_col]] %in% neighborhood, , drop = FALSE]
  }

  ## drop squirrel monkey contamination if present
  if ("species_label" %in% names(anno)) {
    anno <- anno[anno$species_label != "S.sciureus", , drop = FALSE]
  }

  ## ---- counts --------------------------------------------------------------
  if (is.null(counts_file)) {
    counts_file <- if (file.exists(file.path(path, "counts.feather"))) "counts.feather" else "data.feather"
  }
  counts <- feather::read_feather(file.path(path, counts_file))
  counts <- as.data.frame(counts)

  if ("gene" %in% names(counts)) {
    rownames(counts) <- counts$gene
    counts$gene <- NULL
  } else {
    ## Patch-Seq style: sample_id is the row key, transpose to genes x cells
    rownames(counts) <- counts$sample_id
    counts$sample_id <- NULL
    counts <- as.data.frame(t(counts))
  }

  ## align cells to surviving anno rows
  keep <- match(anno$sample_id, colnames(counts), nomatch = 0)
  counts <- counts[, keep[keep > 0], drop = FALSE]
  anno   <- anno[keep > 0, , drop = FALSE]

  if (!identical(colnames(counts), anno$sample_id)) {
    stop("Counts colnames and anno$sample_id do not align after filtering.")
  }

  ## ---- normalize ----------------------------------------------------------
  if (normalize == "cpm_log2") {
    csum <- colSums(counts)
    csum[csum == 0] <- NA
    counts <- sweep(counts, 2, csum, FUN = "/") * 1e6
  }
  counts <- log2(counts + 1)

  ## ---- patch-seq metadata join --------------------------------------------
  meta <- anno
  rownames(meta) <- meta$sample_id

  ## normalize AIBS taxonomy *_label columns (10x + PatchSeq) so subclass /
  ## class / neighborhood / cluster identities match the hier convention
  ## ("L5 ET" → "L5_ET"). Idempotent; safe on cohorts that lack these cols.
  meta <- proj5HT::normalize_taxonomy_labels(meta)

  if ("cell_name_label" %in% names(meta)) {
    md <- tryCatch(headBCI::sheets$cleanMD, error = function(e) NULL)
    if (!is.null(md) && "cell_name" %in% names(md)) {
      meta$cell_name <- meta$cell_name_label
      meta <- dplyr::left_join(meta, md, by = "cell_name", suffix = c("", ".md"))
      ## normalize hier_*_assignment labels in place
      ## (single source of truth = proj5HT::normalize_hier_assignments)
      meta <- proj5HT::normalize_hier_assignments(meta)
      rownames(meta) <- meta$sample_id
    }
    meta$pipe_origin <- ifelse(is.na(match(meta$cell_name_label,
                                           if (!is.null(md)) md$cell_name else character(0),
                                           nomatch = NA)),
                               "on-pipeline", "off-pipeline")
  }

  ## ---- construct ----------------------------------------------------------
  obj <- Seurat::CreateSeuratObject(counts    = counts,
                                    meta.data = meta,
                                    project   = project)
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT*")

  ## (hier_* columns + normalization already attached via the cleanMD join
  ##  above. add_hier_to_seurat() is reserved for objects built outside this
  ##  builder — see R/hier_metadata5HT.R.)

  obj
}

## ----------------------------------------------------------------------------
## I/O helper
## ----------------------------------------------------------------------------

.save_seurat <- function(obj, name) {
  dir.create(.SEU_OUTDIR, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(.SEU_OUTDIR, paste0(name, ".rds"))
  message("Saving ", path)
  saveRDS(obj, file = path)
  invisible(path)
}

## ----------------------------------------------------------------------------
## Per-cohort wrappers
## ----------------------------------------------------------------------------

#' Off-pipeline (proj5HT cohort) macaque Patch-Seq Seurat object.
#' @export
seurat_offpipeline_macaque <- function(save_object = TRUE) {
  obj <- build_seurat_from_feather(
    path      = .SEU_PATHS$macaque_patchseq,
    normalize = "log2_only",
    project   = "offPipeline_Macaque-PatchSeq"
  )
  if (save_object) {
    .save_seurat(obj, "offPipeline-Macaque_PatchSeq")
    .save_seurat(subset(obj, subset = neighborhood_Corr_label %in% .EXC_NEIGHBORHOODS),
                 "offPipeline-Macaque_PatchSeq_excitatory")
    .save_seurat(subset(obj, subset = !neighborhood_Corr_label %in% .EXC_NEIGHBORHOODS),
                 "offPipeline-Macaque_PatchSeq_inhibitory")
  }
  obj
}

#' Off-pipeline human Patch-Seq Seurat object.
#' @export
seurat_offpipeline_human <- function(save_object = TRUE,
                                     subclass = .EXC_SUBCLASS) {
  obj <- build_seurat_from_feather(
    path      = .SEU_PATHS$human_patchseq,
    subclass  = subclass,
    normalize = "log2_only",
    project   = "opHuman-PatchSeq"
  )
  if (save_object) .save_seurat(obj, "offPipeline-Human_PatchSeq")
  obj
}

#' Full macaque Patch-Seq (excitatory) Seurat object.
#' @export
seurat_exc_patchseq_macaque <- function(save_object = TRUE,
                                        subclass = .EXC_SUBCLASS) {
  obj <- build_seurat_from_feather(
    path      = .SEU_PATHS$macaque_patchseq,
    subclass  = subclass,
    normalize = "log2_only",
    project   = "excMacaque-PatchSeq"
  )
  if (save_object) .save_seurat(obj, "excMacaque_PatchSeq")
  obj
}

#' Full human Patch-Seq (excitatory) Seurat object.
#' @export
seurat_exc_patchseq_human <- function(save_object = TRUE,
                                      subclass = .EXC_SUBCLASS) {
  obj <- build_seurat_from_feather(
    path      = .SEU_PATHS$human_patchseq,
    subclass  = subclass,
    normalize = "log2_only",
    project   = "excHuman Patch-seq"
  )
  if (save_object) .save_seurat(obj, "excHuman_PatchSeq")
  obj
}

#' Macaque 10x FACS Seurat object (all cells, optionally subset by
#' subclass / neighborhood).
#' @export
seurat_macaque_10x <- function(save_object = TRUE,
                               subclass = NULL,
                               neighborhood = NULL,
                               name = "MacaqueMTG_10xFACs") {
  obj <- build_seurat_from_feather(
    path         = .SEU_PATHS$macaque_10x,
    subclass     = subclass,
    neighborhood = neighborhood,
    normalize    = "cpm_log2",
    project      = "Macaque10xFACS"
  )
  if (save_object) .save_seurat(obj, name)
  obj
}

#' Human 10x FACS Seurat object.
#' @export
seurat_human_10x <- function(save_object = TRUE,
                             subclass = .EXC_SUBCLASS) {
  obj <- build_seurat_from_feather(
    path        = .SEU_PATHS$human_10x,
    subclass    = subclass,
    normalize   = "cpm_log2",
    project     = "Human10xFACS"
  )
  if (save_object) .save_seurat(obj, "HumanMTG_10xFACs")
  obj
}

#' Chimp 10x FACS Seurat object.
#' @export
seurat_chimp_10x <- function(save_object = TRUE,
                             subclass = c("L5 IT", "L5 ET", "L2/3 IT")) {
  obj <- build_seurat_from_feather(
    path        = .SEU_PATHS$chimp_10x,
    subclass    = subclass,
    normalize   = "cpm_log2",
    project     = "Chimp10xFACS"
  )
  if (save_object) .save_seurat(obj, "ChimpMTG_10xFACs")
  obj
}

#' Gorilla 10x FACS Seurat object.
#' @export
seurat_gorilla_10x <- function(save_object = TRUE,
                               subclass = c("L5 IT", "L5 ET", "L2/3 IT")) {
  obj <- build_seurat_from_feather(
    path        = .SEU_PATHS$gorilla_10x,
    subclass    = subclass,
    normalize   = "cpm_log2",
    project     = "Gorilla10xFACS"
  )
  if (save_object) .save_seurat(obj, "GorillaMTG_10xFACs")
  obj
}

## ----------------------------------------------------------------------------
## Driver
## ----------------------------------------------------------------------------

#' Build every Seurat cohort in one call. Returns a named list.
#' @export
build_all_seurat_objects <- function(save_object = TRUE) {
  list(
    macaque_patchseq_off = seurat_offpipeline_macaque(save_object),
    human_patchseq_off   = seurat_offpipeline_human(save_object),
    macaque_patchseq_exc = seurat_exc_patchseq_macaque(save_object),
    human_patchseq_exc   = seurat_exc_patchseq_human(save_object),
    macaque_10x_all      = seurat_macaque_10x(save_object,
                                              subclass = .EXC_INH_SUBCLASS),
    macaque_10x_exc      = seurat_macaque_10x(save_object,
                                              neighborhood = .EXC_NEIGHBORHOODS,
                                              name = "MacaqueMTG_10xFACs_excitatory"),
    macaque_10x_inh      = seurat_macaque_10x(save_object,
                                              neighborhood = setdiff(
                                                unique_neighborhoods_macaque(),
                                                c(.EXC_NEIGHBORHOODS, "glia")),
                                              name = "MacaqueMTG_10xFACs_inhibitory"),
    human_10x            = seurat_human_10x(save_object),
    chimp_10x            = seurat_chimp_10x(save_object),
    gorilla_10x          = seurat_gorilla_10x(save_object)
  )
}

## helper used by build_all_seurat_objects() for the inhibitory subset
unique_neighborhoods_macaque <- function() {
  a <- feather::read_feather(file.path(.SEU_PATHS$macaque_10x, "anno.feather"))
  unique(a$neighborhood_label)
}

## ----------------------------------------------------------------------------
## (No auto-execution at source-time — call build_all_seurat_objects()
##  explicitly when you want a rebuild.)
## ----------------------------------------------------------------------------
