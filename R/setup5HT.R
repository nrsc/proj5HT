# ============================================================
# setup5HT.R — Project symlink setup for proj5HT
# ============================================================
#
# Thin wrapper around headBCI::bci_setup_links() with
# proj5HT-specific network paths pre-configured.
#
# Usage:
#   proj5HT::setup_links5HT()               # create symlinks
#   proj5HT::setup_links5HT(dry_run = TRUE)  # preview only
# ============================================================

#' Set up proj5HT data symlinks
#'
#' Creates symbolic links from Allen network storage to
#' the local \code{den/} and \code{rookery} directories.
#' Works on Linux, macOS, and Windows.
#'
#' @param dry_run Logical. If TRUE, show what would be created
#'   without actually making links.
#' @param verbose Logical. Print status messages (default TRUE).
#'
#' @return A data.frame of link results (invisibly).
#' @export
setup_links5HT <- function(dry_run = FALSE, verbose = TRUE) {
  if (!requireNamespace("headBCI", quietly = TRUE)) {
    stop("headBCI package is required. Install with: devtools::install('~/headBCI')",
         call. = FALSE)
  }

  proj_root <- path.expand(file.path("~", "proj5HT"))

  links <- headBCI::bci_link_spec(
    den_sources = c(
      "/allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/NHP_expts",
      "/allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/Human_expts",
      "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/10x_seq/GreatApes_Macaque_NCBI",
      "/allen/programs/celltypes/workgroups/hct/HCT_Ephys_Data/manuscript_code/serotonin"
    ),
    den_dest       = file.path(proj_root, "den"),
    rookery_source = "/allen/programs/celltypes/workgroups/hct/SawchukS/rookery",
    rookery_dest   = file.path(proj_root, "rookery")
  )

  if (verbose) message("Setting up proj5HT symlinks...")
  headBCI::bci_setup_links(links = links, dry_run = dry_run, verbose = verbose)
}
