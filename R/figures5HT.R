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
    xlim  = c(-10, 50),
    ylim  = c(0, 200),
    title = "5CT response across pyramids"
  )

  # Ketanserin bath
  figs$ket <- plot_puff_traces(
    dfs$df_ket,
    subclasses = pyr3,
    xlim       = c(-20, 100),
    ylim       = c(0, 200),
    title      = "5HT puff under Ket bath"
  )

  # WAY-100635 bath
  figs$way <- plot_puff_traces(
    dfs$df_way,
    xlim  = c(-20, 80),
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
    xlim       = c(-50, 50),
    ylim       = c(50, 150),
    title      = "Control puff experiments"
  )

  # Ketanserin contamination window
  figs$kcon <- plot_puff_traces(
    dfs$df_kcon,
    xlim  = c(-10, 60),
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
    rapid_window    = c(0, 10),
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
