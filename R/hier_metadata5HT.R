## R/hier_metadata5HT.R
## -----------------------------------------------------------------------------
## Hierarchical-classifier (`hier`) metadata helpers, sourced from
## `headBCI::sheets$cleanMD`.
##
## The hier_* columns are built by the headBCI sheets pipeline (28 columns,
## covering class / neighborhood / subclass / cluster + their bootstrap /
## aggregate / avg-correlation / runner-up / directly-assigned variants). This
## file exposes:
##
##   hier5HT()                 per-cell hier metadata (cell_name keyed),
##                             with `hier_subclass_assignment` normalized to
##                             the proj5HT cohort convention (e.g. "L2/3 IT"
##                             → "L2/3_IT").
##
##   add_hier_to_seurat(obj)   left-joins hier5HT() onto obj@meta.data via
##                             cell_name. Adds the hier_* columns; existing
##                             columns are left untouched.
## -----------------------------------------------------------------------------

#' Per-cell hierarchical-classifier metadata.
#'
#' Pulls hier_* columns from `headBCI::sheets$cleanMD`, normalizes
#' subclass/cluster strings to the proj5HT convention (spaces replaced with
#' underscores), and returns one row per `cell_name`.
#'
#' @param cols character vector of hier_* columns to return. Defaults to the
#'   most commonly used subset; pass `NULL` to return all hier_* columns.
#' @return data.frame with at least `cell_name` plus the requested hier columns.
#' Normalize categorical label columns to the proj5HT convention.
#'
#' Core string normalizer: replaces spaces with underscores
#' (e.g. `"L2/3 IT"` → `"L2/3_IT"`) on any column whose name is in `cols`
#' and is present in `df`. Empty strings and `NA`s pass through unchanged.
#' Idempotent.
#'
#' This is the single source of truth used by both `normalize_hier_assignments()`
#' (headBCI hier classifier output) and `normalize_taxonomy_labels()` (AIBS
#' 10x / PatchSeq taxonomy `*_label` columns).
#'
#' @param df data.frame (or any 2D structure accepting `df[[col]] <- ...`).
#' @param cols character vector of column names to normalize.
#' @return `df` with the named columns normalized in place.
#' @export
normalize_label_cols <- function(df, cols) {
  for (c in intersect(cols, colnames(df))) {
    v <- df[[c]]
    df[[c]] <- ifelse(is.na(v) | !nzchar(v), NA_character_,
                      gsub(" ", "_", v, fixed = TRUE))
  }
  df
}

#' Normalize `hier_*_assignment` labels (headBCI hier classifier).
#'
#' Thin wrapper over `normalize_label_cols()` with the four canonical hier
#' assignment columns as the default. Use anywhere hier assignment labels are
#' joined in from `headBCI::sheets$cleanMD`.
#'
#' @param df data.frame.
#' @param cols character vector of column names. Defaults to the four
#'   `hier_*_assignment` columns.
#' @return `df` with the named columns normalized in place.
#' @export
normalize_hier_assignments <- function(df,
                                       cols = c("hier_class_assignment",
                                                "hier_neighborhood_assignment",
                                                "hier_subclass_assignment",
                                                "hier_cluster_assignment")) {
  normalize_label_cols(df, cols)
}

#' Normalize AIBS taxonomy `*_label` columns (10x reference + PatchSeq).
#'
#' Allen Brain Atlas 10x feathers and PatchSeq annotation feathers carry
#' subclass / class / neighborhood / cluster identities under `*_label`
#' columns (and `*_Corr_label` variants on PatchSeq) with space-separated
#' values (`"L5 ET"`, `"L2/3 IT"`). Normalize them to the underscored form
#' so axis text and equality filters match the hier convention across
#' panels.
#'
#' @param df data.frame.
#' @param cols character vector of column names. Defaults to the canonical
#'   AIBS taxonomy label columns plus their `_Corr_label` variants.
#' @return `df` with the named columns normalized in place.
#' @export
normalize_taxonomy_labels <- function(df,
                                      cols = c("class_label",
                                               "neighborhood_label",
                                               "subclass_label",
                                               "cluster_label",
                                               "class_Corr_label",
                                               "neighborhood_Corr_label",
                                               "subclass_Corr_label",
                                               "cluster_Corr_label")) {
  normalize_label_cols(df, cols)
}

#' @export
hier5HT <- function(cols = c("hier_class_assignment",
                             "hier_neighborhood_assignment",
                             "hier_subclass_assignment",
                             "hier_cluster_assignment",
                             "hier_subclass_bootstrap_prob",
                             "hier_subclass_avg_corr",
                             "hier_class_bootstrap_prob",
                             "hier_cluster_bootstrap_prob")) {
  cmd <- headBCI::sheets$cleanMD
  if (!"cell_name" %in% colnames(cmd)) {
    stop("headBCI::sheets$cleanMD has no cell_name column")
  }

  hier_cols_all <- grep("^hier_", colnames(cmd), value = TRUE)
  if (is.null(cols)) cols <- hier_cols_all
  missing <- setdiff(cols, hier_cols_all)
  if (length(missing)) {
    warning("hier columns missing from headBCI::sheets$cleanMD (filling NA): ",
            paste(missing, collapse = ", "))
    for (c in missing) cmd[[c]] <- NA
  }

  out <- cmd[, c("cell_name", cols), drop = FALSE]
  out <- normalize_hier_assignments(out)

  ## one row per cell (cleanMD is already cell-level, but guard anyway)
  out <- out[!is.na(out$cell_name) & nzchar(out$cell_name), ]
  out <- out[!duplicated(out$cell_name), ]
  rownames(out) <- NULL
  out
}

#' Overlay hierarchical-classifier metadata onto a Seurat object.
#'
#' Convenience wrapper for objects not built via
#' `dev/Seurat/Seurat_Objects.R::build_seurat_from_feather()` (which already
#' joins `cleanMD` and normalizes in-place). Left-joins `hier5HT()` to
#' `obj@meta.data` on `cell_name`; any hier columns already present are
#' re-normalized via `normalize_hier_assignments()` instead of being
#' duplicated.
#'
#' @param obj Seurat object whose `@meta.data` contains `cell_name`.
#' @param cols character vector of hier columns to attach (see `hier5HT()`).
#' @return the same object with hier columns appended / normalized in
#'   `@meta.data`.
#' @export
add_hier_to_seurat <- function(obj, cols = NULL) {
  if (!"cell_name" %in% colnames(obj@meta.data)) {
    stop("obj@meta.data has no `cell_name` column — cannot join hier metadata")
  }
  hier <- if (is.null(cols)) hier5HT(cols = NULL) else hier5HT(cols = cols)
  hier_cols <- setdiff(colnames(hier), "cell_name")

  md <- obj@meta.data

  ## drop any pre-existing hier columns so the freshly-normalized values win
  pre_existing <- intersect(hier_cols, colnames(md))
  if (length(pre_existing)) md[pre_existing] <- NULL

  md$.rowkey <- rownames(md)
  md <- merge(md, hier, by = "cell_name", all.x = TRUE, sort = FALSE)
  rownames(md) <- md$.rowkey
  md$.rowkey <- NULL
  md <- md[colnames(obj), , drop = FALSE]   ## restore Seurat row ordering
  obj@meta.data <- md

  matched <- sum(!is.na(obj@meta.data$hier_subclass_assignment))
  message(sprintf("hier metadata attached: %d / %d cells matched",
                  matched, ncol(obj)))
  obj
}
