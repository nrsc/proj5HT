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
