# ============================================================
# params5HT.R — Project-level parameters for proj5HT
# ============================================================
#
# Follows the same pattern as projHCT::paramsHCT().
#
# The `params5HT` list is saved to `data/params5HT.rda` so it is
# available as `proj5HT::params5HT` via LazyData.  Call
# `paramsHCT5HT()` to regenerate and re-install after edits.
#
# Usage in other functions:
#   mMD_path default → params5HT$paths$master_ui
#   mMD_sheet default → params5HT$sheets$master_macaque
# ============================================================

#' Rebuild and install proj5HT project parameters
#'
#' Constructs the canonical \code{params5HT} list, compares it to the
#' installed version, and (if different) updates \code{data/params5HT.rda},
#' re-documents, and re-installs the package.
#'
#' @param return logical; if TRUE, silently return the params list.
#'
#' @return The \code{params5HT} list (invisibly unless \code{return = TRUE}).
#' @export
params5HT_build <- function(return = FALSE) {

  params5HT <- list(

    # ---- project identity ----
    home = "~",
    proj = "proj5HT",

    # ---- paths ----
    paths = list(
      # Canonical root for the project (tilde-safe)
      project_root  = file.path("~", "proj5HT"),
      rookery       = file.path("~", "proj5HT", "rookery"),

      # ODS metadata sheets
      master_ui     = file.path("~", "proj5HT", "den", "serotonin", "Data", "Master_UI_data.ods"),
      selected_md   = file.path("~", "proj5HT", "den", "serotonin", "Data", "selected_merged_metadata.ods"),

      # Compiled data
      compiled_rds  = file.path("~", "proj5HT", "data-raw", "compiled5HT.rds"),
      compile_puff  = file.path("~", "proj5HT", "data-raw", "compile_srtPuff.csv"),
      compile_dwash = file.path("~", "proj5HT", "data-raw", "compile_srtDrugWash.csv"),
      compile_kwash = file.path("~", "proj5HT", "data-raw", "compile_srtKetWash.csv"),

      # Data sources (Allen paths — will be symlinked via Makefile)
      den_NHP    = "//allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/NHP_expts",
      den_Human  = "//allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/Human_expts",
      den_srt    = "//allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/manuscript_code/serotonin"
    ),

    # ---- ODS sheet names ----
    sheets = list(
      master_macaque  = "Master_Macaque",
      master_human    = "Human",
      selected_mac    = "Macaque",
      selected_human  = "Human"
    ),

    # ---- nnest tags ----
    tags = list(
      parent  = "-bci.rds",
      child   = "-srt.rds",
      puff    = "-srtPuff.csv",
      dwash   = "-srtDrugWash.csv",
      kwash   = "-srtKetWash.csv",
      rmp     = "-rmp.csv",
      asp     = "-asp.csv",
      md_snap = "-5ht_md_snapshot.csv"
    ),

    # ---- classifier defaults ----
    classifier = list(
      baseline_window   = c(-5, 0),
      rapid_window      = c(0, 10),
      mid_window        = c(10, 40),
      late_window       = c(40, 60),
      post_window       = c(0, 60),
      min_excursion_pct = 10,
      min_persist_s     = 0.5,
      smooth_method     = "rollmean",
      smooth_k          = 3,
      loess_span        = 0.2,
      subclass_keep     = c("L23_IT", "L5_ET", "L5_IT"),
      recode_subclass   = c("L3c" = "L23_IT"),
      bath_keep         = "none",
      puff_pattern      = "^5HT"
    ),

    # ---- ASP / build defaults ----
    asp = list(
      ttl_default = 5,
      outlier_sd  = 5,
      bl_window   = c(-5, 0)
    ),

    # ---- cells excluded from analysis ----
    exclude = c(
      "Q21.26.010.11.13.03",
      "Q21.26.014.1A.21.02",
      "Q21.26.010.11.13.02"
    ),

    # ---- version tracking ----
    version = list(
      params_version = "0.1.0"
    )
  )

  # ---- check whether installed copy is current ----
  if (requireNamespace("proj5HT", quietly = TRUE) &&
      exists("params5HT", where = asNamespace("proj5HT"), inherits = FALSE)) {
    installed <- get("params5HT", envir = asNamespace("proj5HT"))
    if (identical(installed, params5HT)) {
      message("params5HT is already up to date.")
      if (return) return(params5HT)
      return(invisible(params5HT))
    }
  }

  # ---- save, document, and re-install ----
  message("Updating params5HT dataset for proj5HT package")
  usethis::use_data(params5HT, overwrite = TRUE)
  devtools::document()
  system(
    paste(
      "R CMD INSTALL --preclean --no-multiarch --with-keep.source",
      path.expand(file.path(params5HT$home, params5HT$proj))
    )
  )

  if (return) return(params5HT)
  invisible(params5HT)
}


#' @title proj5HT project parameters
#'
#' @description
#' A list of project-level parameters: canonical paths, ODS sheet names,
#' nnest tags, classifier defaults, ASP defaults, and excluded cells.
#'
#' Access as \code{proj5HT::params5HT} (lazy-loaded from
#' \code{data/params5HT.rda}).  Rebuild with \code{params5HT_build()}.
#'
#' @format A named list.
#' @examples
#' proj5HT::params5HT$paths$master_ui
#' proj5HT::params5HT$classifier$rapid_window
"params5HT"
