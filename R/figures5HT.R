#' proj5HT figure orchestrator
#'
#' Builds the standard family of dataframes via `figure_df_build()` and
#' renders the project's stock figures using `plot_puff_traces()` for a
#' consistent style (faint grey per-cell traces, black averaged trace on top).
#'
#' The function returns the list of `ggplot` objects so callers can print,
#' arrange or save them. Set `interactive = TRUE` to also open the underlying
#' figure-script files for interactive editing.
#'
#' @param dfs          Optional pre-built list from `figure_df_build()`.
#' @param interactive  If TRUE, opens auxiliary figure scripts via `file.edit`.
#'
#' @return Named list of `ggplot` objects.
#' @export
figures5HT <- function(dfs = NULL, interactive = FALSE) {

  ## -------------------------------------------------------------------- setup
  source("~/proj5HT/R/figures_df_build.R")
  source("~/proj5HT/figs/scripts/plot_helpers.R")
  if (is.null(dfs)) dfs <- figure_df_build()

  pyr3 <- c("L23_IT", "L5_IT", "L5_ET")
  pyr4 <- c("L23_IT", "L5_IT", "L5_ET", "L4_IT")

  ## ----------------------------------------------------------------- figures
  figs <- list()

  # Standard 5HT puff response across macaque pyramids (df0)
  figs$std_puff <- plot_puff_traces(
    dfs$df0,
    subclasses = pyr3,
    xlim       = c(-10, 50),
    ylim       = c(0, 300),
    title      = "5HT puff response — macaque pyramids"
  )

  # L23_IT split by cortical depth (superficial < 500 µm <= deep)
  figs$l23_depth <- plot_puff_traces(
    dplyr::filter(dfs$df0,
                  assigned_subclass == "L23_IT",
                  !is.na(depth_bin)),
    facet = "depth_bin",
    xlim  = c(-10, 50),
    ylim  = c(0, 300),
    title = "L23_IT — superficial vs deep (split at 500 µm)"
  )

  # L23_IT split by transcriptomic cluster_Corr (Exc L2/3 clusters), with
  # per-cell traces coloured by cortical depth (superficial vs deep).
  l23_clusters <- c("Exc L2-3 LINC00507 FREM3",
                    "Exc L2 LAMP5 LTK",
                    "L2/3_IT_1",
                    "L2/3_IT_3")
  depth_palette <- c(superficial = "#1f77b4", deep = "#d62728")
  figs$l23_cluster <- plot_puff_traces(
    dplyr::filter(dfs$df0,
                  assigned_subclass == "L23_IT",
                  cluster_Corr %in% l23_clusters,
                  !is.na(depth_bin)),
    facet    = "cluster_Corr",
    color_by = "depth_bin",
    palette  = depth_palette,
    xlim     = c(-10, 50),
    ylim     = c(0, 300),
    title    = "L23_IT — split by cluster_Corr, coloured by depth"
  )

  # Human 5HT puff response
  figs$human <- plot_puff_traces(
    dfs$dfH,
    subclasses = pyr4,
    xlim       = c(-10, 50),
    ylim       = c(0, 300),
    title      = "5HT response across human cortical pyramids"
  )

  # 5CT puff
  figs$puff_5ct <- plot_puff_traces(
    dfs$df_5ct,
    subclasses = pyr3,
    xlim  = c(-10, 50),
    ylim  = c(0, 200),
    title = "5CT response across pyramids"
  )

  # Ketanserin bath
  figs$ket <- plot_puff_traces(
    dfs$df_ket,
    subclasses = pyr3,
    xlim       = c(-10, 50),
    ylim       = c(0, 200),
    title      = "5HT puff under Ket bath"
  )

  # WAY-100635 bath
  figs$way <- plot_puff_traces(
    dfs$df_way,
    xlim  = c(-10, 50),
    ylim  = c(0, 300),
    title = "5HT/5CT response under WAY-100635"
  )

  # Carbachol
  figs$carb <- plot_puff_traces(
    dfs$df_carb,
    xlim  = c(-10, 60),
    ylim  = c(0, 300),
    title = "Carbachol puff response"
  )

  # Control (acsf / control puff) — sanity check
  figs$control <- plot_puff_traces(
    dfs$df_control,
    subclasses = pyr3,
    xlim       = c(-20, 50),
    ylim       = c(50, 150),
    title      = "Control puff experiments"
  )

  # Ketanserin contamination window
  figs$kcon <- plot_puff_traces(
    dfs$df_kcon,
    xlim  = c(-10, 50),
    ylim  = c(0, 300),
    title = "Cells run during Ket contamination window"
  )

  # Inhibitory subclasses
  figs$inh <- plot_puff_traces(
    dfs$df_inh,
    xlim  = c(-10, 50),
    ylim  = c(0, 300),
    title = "5HT response — inhibitory subclasses"
  )

  ## ----------------------------------------------- 4-way response classifier
  per_cell <- classify_response5HT(
    dplyr::filter(dfs$df0, assigned_subclass %in% pyr3),
    rapid_window    = c(-1, 10),
    extended_window = c(10, 50),
    thr_pct         = 15,
    min_persist_s   = 1,
    smooth_k        = 5
  )
  figs$per_cell <- per_cell

  figs$by_response <- plot_puff_traces_by_response(
    dfs$df0, per_cell,
    xlim  = c(-10, 50),
    ylim  = c(0, 300),
    title = "5HT puff — grouped by response classification"
  )

  # Per-cell inspector for QC of the "excitation" panel
  figs$exc_inspector <- plot_excitatory_inspector(
    dfs$df0, per_cell,
    n_panels = 4,
    xlim     = c(-10, 50),
    ylim     = c(0, 350),
    title    = "Classified excitation — per-cell QC"
  )

  # Diagnostic table: deepest dip inside each excitation-classified cell.
  # Use this to spot cells whose silent period sits outside the classifier's
  # rapid window (e.g. a 5 s metadata offset would put the dip at t ~ 5).
  figs$exc_dips <- diagnose_excitation_dips(
    dfs$df0, per_cell,
    window       = c(-2, 15),
    zero_thr_pct = 5
  )

  ## ------------------------------------------------- auxiliary figure scripts
  if (interactive) {
    file.edit("~/proj5HT/R/figures_df_build.R")
    file.edit("~/proj5HT/figs/scripts/plot_helpers.R")
    file.edit("~/proj5HT/figs/scripts/srtHumanData.R")
    file.edit("~/proj5HT/figs/scripts/srtFigDateCompare.R")
    file.edit("~/proj5HT/figs/scripts/srt_ketContam.R")
    file.edit("~/proj5HT/figs/scripts/srtFigKet.R")
    file.edit("~/proj5HT/figs/scripts/srtFig5CT.R")
    file.edit("~/proj5HT/figs/scripts/srtFig5CTpuff.R")
  }

  invisible(figs)
}
