#' Styled puff trace figure: grey per-cell traces + black average
#'
#' Standard figure style used across proj5HT puff figures. Each cell is drawn
#' as a faint grey line; the across-cell mean (per `x_bin` per facet) is
#' overlaid in black. Facets default to `assigned_subclass`.
#'
#' Expected columns on `df`: `x_bin`, `cell_name`, `percent_change`, and the
#' facet variable (default `assigned_subclass`).
#'
#' @param df         Long-format dataframe (e.g. one of the dfs from
#'                   `figure_df_build()`).
#' @param subclasses Optional character vector of subclasses to keep.
#' @param facet      Faceting variable name (string). Default
#'                   `"assigned_subclass"`. Set to `NA` to disable faceting.
#' @param xlim,ylim  Axis limits passed to `coord_cartesian()`.
#' @param title      Plot title.
#' @param avg_fun    Summary function for the black overlay (default `mean`).
#' @param trace_alpha,trace_lwd   Styling for individual cell traces.
#' @param mean_lwd               Line width for the black average.
#' @param min_cells_per_bin      Bins with fewer cells than this are dropped
#'                               from the average overlay.
#'
#' @return A `ggplot` object.
#' @export
plot_puff_traces <- function(df,
                             subclasses        = NULL,
                             facet             = "assigned_subclass",
                             xlim              = c(-10, 50),
                             ylim              = c(0, 300),
                             title             = NULL,
                             avg_fun           = mean,
                             trace_alpha       = 0.35,
                             trace_lwd         = 0.4,
                             mean_lwd          = 1.1,
                             bin_step          = 2,
                             min_cells_per_bin = 2) {

  if (!is.null(subclasses)) {
    df <- dplyr::filter(df, assigned_subclass %in% subclasses)
  }

  # Restrict to the plot window so completion below uses a sensible grid.
  df <- dplyr::filter(df, x_bin >= xlim[1], x_bin <= xlim[2])

  group_cols <- c("x_bin", "cell_name")
  if (!is.na(facet)) group_cols <- c(group_cols, facet)

  # Per-cell, per-bin mean of `percent_change`.
  per_cell <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(y = avg_fun(percent_change, na.rm = TRUE),
                     .groups = "drop")

  # Ensure each cell has every bin in the plot window — bins with no spikes
  # get y = 0 so they contribute to the across-cell average rather than being
  # silently dropped.
  grid <- seq(xlim[1], xlim[2], by = bin_step)
  cell_keys <- if (!is.na(facet)) c("cell_name", facet) else "cell_name"
  per_cell <- per_cell %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(cell_keys))) %>%
    tidyr::complete(x_bin = grid, fill = list(y = 0)) %>%
    dplyr::ungroup()

  mean_group_cols <- c("x_bin")
  if (!is.na(facet)) mean_group_cols <- c(mean_group_cols, facet)

  df_mean <- per_cell %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(mean_group_cols))) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      y       = avg_fun(y, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cells >= min_cells_per_bin)

  p <- ggplot2::ggplot(per_cell,
                       ggplot2::aes(x = x_bin, y = y, group = cell_name)) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2,
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_line(colour = "grey60",
                       alpha    = trace_alpha,
                       linewidth = trace_lwd,
                       na.rm    = TRUE) +
    ggplot2::geom_line(data        = df_mean,
                       mapping     = ggplot2::aes(x = x_bin, y = y),
                       inherit.aes = FALSE,
                       colour      = "black",
                       linewidth   = mean_lwd,
                       na.rm       = TRUE) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(title = title,
                  x     = "Time (s)",
                  y     = "Percent change") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   panel.grid.minor = ggplot2::element_blank())

  if (!is.na(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  }

  p
}


#' Plot puff traces faceted by 4-way response classification
#'
#' Convenience wrapper that joins the per-cell `response_type` from
#' `classify_response5HT()` back onto a long-format trace dataframe and
#' renders the standard grey-traces / black-mean style faceted by class.
#'
#' @param df       Long-format dataframe (one of the dfs from
#'                 `figure_df_build()`).
#' @param per_cell tibble from `classify_response5HT()`.
#' @param drop_no_change  drop the `no_change` panel (default `TRUE`).
#' @param ...      passed through to `plot_puff_traces()`.
#' @return ggplot
#' @export
plot_puff_traces_by_response <- function(df, per_cell,
                                         drop_no_change = TRUE,
                                         title = "Responses by class",
                                         ...) {
  df2 <- dplyr::left_join(
    df,
    dplyr::select(per_cell, cell_name, response_type),
    by = "cell_name"
  ) %>%
    dplyr::filter(!is.na(response_type))

  if (drop_no_change) {
    df2 <- dplyr::filter(df2, response_type != "no_change")
  }

  plot_puff_traces(df2,
                   facet = "response_type",
                   title = title,
                   ...)
}


#' Inspector plot: classified-excitatory cells across N stacked panels
#'
#' Builds a stacked, vertically-faceted plot of every cell whose classifier
#' call is `"excitation"`. Cells are split into `n_panels` groups so each
#' panel holds a manageable number of traces, every cell gets a distinct
#' colour, and the legend lets you identify offending traces by `cell_name`.
#'
#' Intended for QC/spot-checking the classifier — not a publication figure.
#'
#' @param df         Long-format trace dataframe (e.g. `dfs$df0`).
#' @param per_cell   Output of `classify_response5HT()`.
#' @param n_panels   Number of stacked panels (default `4`).
#' @param xlim,ylim  Axis limits.
#' @param title      Plot title.
#' @param bin_step   Bin width used to summarise per cell (default `2` s).
#'
#' @return ggplot object.
#' @export
plot_excitatory_inspector <- function(df, per_cell,
                                      n_panels = 4,
                                      xlim     = c(-10, 50),
                                      ylim     = c(0, 350),
                                      title    = "Classified excitation — per-cell inspector",
                                      bin_step = 2) {

  exc_cells <- per_cell %>%
    dplyr::filter(as.character(response_type) == "excitation") %>%
    dplyr::pull(cell_name) %>%
    unique() %>%
    sort()

  if (length(exc_cells) == 0) {
    stop("No cells classified as 'excitation' in per_cell.")
  }

  # Even split into n_panels groups (alphabetical for reproducibility).
  panel_assignment <- tibble::tibble(
    cell_name = exc_cells,
    panel     = paste0(
      "panel ",
      as.integer(cut(seq_along(exc_cells),
                     breaks = n_panels,
                     labels = FALSE,
                     include.lowest = TRUE))
    )
  )

  d <- df %>%
    dplyr::filter(cell_name %in% exc_cells,
                  x_bin >= xlim[1], x_bin <= xlim[2]) %>%
    dplyr::group_by(cell_name, x_bin) %>%
    dplyr::summarise(y = mean(percent_change, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::left_join(panel_assignment, by = "cell_name")

  # Per-panel palette so adjacent cells in the same panel are distinguishable.
  d <- d %>% dplyr::mutate(cell_name = factor(cell_name, levels = exc_cells))

  ggplot2::ggplot(d, ggplot2::aes(x = x_bin, y = y,
                                  group = cell_name,
                                  colour = cell_name)) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2,
                        colour = "grey60", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.7, na.rm = TRUE) +
    ggplot2::facet_wrap(~ panel, ncol = 1) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::scale_colour_viridis_d(option = "turbo") +
    ggplot2::labs(title = title,
                  x      = "Time (s)",
                  y      = "Percent change",
                  colour = "cell_name") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      legend.text      = ggplot2::element_text(size = 7),
      legend.key.height = ggplot2::unit(0.7, "lines")
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1))
}


#' Diagnose deep dips inside "excitation"-classified cells
#'
#' Returns a tibble (one row per cell currently classified as
#' `"excitation"`) reporting the minimum raw `percent_change` and the
#' time at which it occurred, looking inside a configurable window. This
#' is the fastest way to see whether an apparent silent period in the
#' classifier panel is sitting in the expected rapid window (t ~ 0) or
#' whether it has been shifted (e.g. a 5 s metadata offset would put the
#' minimum at t ~ 5 rather than t ~ 0).
#'
#' Cells with `min_pct <= zero_thr_pct` *should* have been flagged
#' biphasic — if they slipped through, the window almost certainly does
#' not cover where their silence actually lies.
#'
#' @param df             Long-format trace dataframe (e.g. `dfs$df0`).
#' @param per_cell       Output of `classify_response5HT()`.
#' @param window         Time range to scan for the dip (default
#'                       `c(-2, 15)` — wider than the classifier's rapid
#'                       window on purpose, to catch shifted onsets).
#' @param zero_thr_pct   Threshold used to flag a "deep" dip (default `5`).
#' @param sort_by_min    If `TRUE` (default) rows are sorted by ascending
#'                       `min_pct` so the worst offenders are at the top.
#'
#' @return tibble with columns: `cell_name`, `assigned_subclass`,
#'   `min_pct`, `min_time`, `is_deep_dip`, plus any per-cell metadata
#'   carried by `per_cell`.
#' @export
diagnose_excitation_dips <- function(df, per_cell,
                                     window       = c(-2, 15),
                                     zero_thr_pct = 5,
                                     sort_by_min  = TRUE) {

  exc_cells <- per_cell %>%
    dplyr::filter(as.character(response_type) == "excitation") %>%
    dplyr::pull(cell_name) %>%
    unique()

  if (length(exc_cells) == 0) {
    return(tibble::tibble())
  }

  out <- df %>%
    dplyr::filter(cell_name %in% exc_cells,
                  time >= window[1], time <= window[2]) %>%
    dplyr::group_by(cell_name) %>%
    dplyr::summarise(
      min_pct  = suppressWarnings(min(percent_change, na.rm = TRUE)),
      min_time = time[which.min(percent_change)],
      .groups  = "drop"
    ) %>%
    dplyr::mutate(is_deep_dip = is.finite(min_pct) & min_pct <= zero_thr_pct)

  # carry useful per-cell context
  carry <- intersect(c("assigned_subclass", "Species", "puff", "bath",
                       "expCon", "blockers", "assigned_depth",
                       "rapid_zero", "ext_zero",
                       "rapid_exc", "rapid_inh", "ext_exc", "ext_inh"),
                     names(per_cell))
  out <- dplyr::left_join(
    out,
    dplyr::select(per_cell, cell_name, dplyr::all_of(carry)),
    by = "cell_name"
  )

  if (sort_by_min) {
    out <- dplyr::arrange(out, min_pct)
  }
  out
}
