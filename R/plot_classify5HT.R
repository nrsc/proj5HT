# ============================================================
# plot_classify5HT.R — Plotting functions for 5-HT response classifier
# ============================================================
#
# Functions here consume the output of classify5HT.R
# (df_prepped + per_cell) and produce ggplot2 figures.
#
# ============================================================


#' QC trace plot for a single cell
#'
#' Plots the raw and smoothed pct_of_baseline trace with shaded
#' time-window bands and excursion thresholds.
#'
#' @param df_prepped data.frame from \code{prep_puff_df()}.
#' @param cell_id character cell_name to plot.
#' @param smooth_method,smooth_k,loess_span smoothing parameters.
#' @param rapid_window,mid_window,late_window time windows shown as shaded regions.
#' @param min_excursion_pct threshold lines drawn at 100 +/- this value.
#'
#' @return ggplot object.
#' @export
plot_cell_qc <- function(df_prepped,
                         cell_id,
                         smooth_method = "rollmean",
                         smooth_k = 3,
                         loess_span = 0.2,
                         rapid_window = c(0, 10),
                         mid_window   = c(10, 40),
                         late_window  = c(40, 60),
                         min_excursion_pct = 10) {

  d <- df_prepped %>%
    dplyr::filter(.data$cell_name == cell_id) %>%
    dplyr::arrange(.data$time)

  if (nrow(d) == 0) stop("No rows for cell_id: ", cell_id)

  d <- d %>%
    dplyr::mutate(
      pct_smooth = smooth_vec(
        y = .data$pct_of_baseline,
        method = smooth_method,
        k = smooth_k,
        loess_span = loess_span,
        x = .data$time
      )
    )

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$time)) +
    ggplot2::annotate("rect",
      xmin = rapid_window[1], xmax = rapid_window[2],
      ymin = -Inf, ymax = Inf, alpha = 0.06
    ) +
    ggplot2::annotate("rect",
      xmin = mid_window[1], xmax = mid_window[2],
      ymin = -Inf, ymax = Inf, alpha = 0.04
    ) +
    ggplot2::annotate("rect",
      xmin = late_window[1], xmax = late_window[2],
      ymin = -Inf, ymax = Inf, alpha = 0.05
    ) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2,
                        linewidth = 0.4, color = "grey40") +
    ggplot2::geom_hline(yintercept = 100 + min_excursion_pct, linetype = 3,
                        linewidth = 0.35, color = "#2166AC") +
    ggplot2::geom_hline(yintercept = 100 - min_excursion_pct, linetype = 3,
                        linewidth = 0.35, color = "#B2182B") +
    ggplot2::geom_line(ggplot2::aes(y = .data$pct_of_baseline),
                       color = "grey70", linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$pct_smooth),
                       color = "black", linewidth = 0.75) +
    ggplot2::labs(
      title    = cell_id,
      subtitle = first_non_missing(d$assigned_subclass),
      x = "Time (s)",
      y = "% baseline"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}


#' Batch QC plots — save many per-cell traces to a directory
#'
#' @param df_prepped data.frame from \code{prep_puff_df()}.
#' @param cell_ids character vector of cell_names to plot.
#' @param out_dir output directory for PNGs.
#' @param width,height,dpi passed to \code{ggplot2::ggsave()}.
#' @param ... extra args forwarded to \code{plot_cell_qc()}.
#'
#' @return invisible NULL. Side effect: writes PNG files.
#' @export
plot_cell_qc_many <- function(df_prepped,
                              cell_ids,
                              out_dir = "qc_bucket_plots",
                              width = 7, height = 4.5, dpi = 200,
                              ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (cid in cell_ids) {
    p  <- plot_cell_qc(df_prepped, cid, ...)
    fn <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", cid), ".png"))
    ggplot2::ggsave(fn, p, width = width, height = height, dpi = dpi)
  }

  message("Saved ", length(cell_ids), " QC plots to: ", out_dir)
  invisible(NULL)
}


#' Stacked bar:  response composition by subclass
#'
#' @param per_cell one-row-per-cell tibble from \code{build_per_cell_classifier()}.
#' @param x column name to group on (default "assigned_subclass").
#'
#' @return ggplot object.
#' @export
plot_response_composition <- function(per_cell, x = "assigned_subclass") {

  assert_has_cols(per_cell, c(x, "response_call"))

  per_cell %>%
    dplyr::filter(!is.na(.data$response_call)) %>%
    dplyr::count(.data[[x]], .data$response_call, name = "n") %>%
    dplyr::group_by(.data[[x]]) %>%
    dplyr::mutate(frac = .data$n / sum(.data$n)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data[[x]], y = .data$frac, fill = .data$response_call
    )) +
    ggplot2::geom_col(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = response_cols, labels = pretty_labels, drop = FALSE
    ) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(x = NULL, y = "Fraction of cells", fill = "Response class") +
    ggplot2::theme_minimal(base_size = 12)
}


#' Tile heatmap:  response call vs subclass
#'
#' @param per_cell one-row-per-cell tibble.
#' @param x,y column names for tile axes.
#'
#' @return ggplot object.
#' @export
plot_response_heatmap <- function(per_cell,
                                  x = "response_call",
                                  y = "assigned_subclass") {

  assert_has_cols(per_cell, c(x, y))

  per_cell %>%
    dplyr::filter(!is.na(.data[[x]]), !is.na(.data[[y]])) %>%
    dplyr::count(.data[[y]], .data[[x]], name = "n") %>%
    dplyr::group_by(.data[[y]]) %>%
    dplyr::mutate(frac = .data$n / sum(.data$n)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(
      x = .data[[x]], y = .data[[y]], fill = .data$frac
    )) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(
      label = scales::percent(.data$frac, accuracy = 1)
    ), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Fraction") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}


#' Violin + jitter of response peak by subclass
#'
#' @param per_cell one-row-per-cell tibble.
#'
#' @return ggplot object.
#' @export
plot_centered_violin <- function(per_cell) {
  assert_has_cols(per_cell, c("assigned_subclass", "response_peak", "direction"))

  ggplot2::ggplot(
    per_cell %>% dplyr::filter(is.finite(.data$response_peak)),
    ggplot2::aes(x = .data$assigned_subclass, y = .data$response_peak)
  ) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_violin(fill = "grey85", color = "black",
                         linewidth = 0.35, trim = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$direction),
      position = ggplot2::position_jitter(width = 0.12),
      size = 1.2, alpha = 0.55
    ) +
    ggplot2::scale_color_manual(
      values = c(
        inhibition = "#B2182B",
        excitation = "#2166AC",
        biphasic   = "purple",
        no_change  = "grey50"
      ),
      guide = "none"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Peak % baseline") +
    ggplot2::theme_minimal(base_size = 12)
}


#' Faceted exc / inh traces with per-cell grey lines and group mean
#'
#' @param df_prepped data.frame from \code{prep_puff_df()}.
#' @param per_cell one-row-per-cell tibble.
#' @param subclasses character vector of subclass labels to include.
#' @param xlim numeric length-2. Time range for x-axis.
#' @param smooth_k window width for per-cell rolling mean.
#'
#' @return ggplot object.
#' @export
plot_exc_inh_traces <- function(df_prepped,
                                per_cell,
                                subclasses = c("L23_IT", "L5_ET", "L5_IT"),
                                xlim = c(-5, 60),
                                smooth_k = 3) {

  assert_has_cols(df_prepped, c("cell_name", "time", "pct_of_baseline"))
  assert_has_cols(per_cell, c("cell_name", "assigned_subclass", "response_type"))

  d0 <- df_prepped %>%
    dplyr::select(-dplyr::any_of(c("assigned_subclass", "response_type"))) %>%
    dplyr::inner_join(
      per_cell %>%
        dplyr::select("cell_name", "assigned_subclass", "response_type"),
      by = "cell_name"
    ) %>%
    dplyr::filter(
      .data$assigned_subclass %in% subclasses,
      .data$response_type %in% c("excitation", "inhibition"),
      .data$time >= xlim[1], .data$time <= xlim[2]
    ) %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(.data$pct_of_baseline,
                                 k = smooth_k, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  d_mean <- d0 %>%
    dplyr::group_by(.data$response_type, .data$assigned_subclass, .data$time) %>%
    dplyr::summarise(
      mean_pct = mean(.data$pct_smooth, na.rm = TRUE),
      .groups = "drop"
    )

  ggplot2::ggplot() +
    ggplot2::facet_grid(response_type ~ assigned_subclass, scales = "free_y") +
    ggplot2::geom_hline(yintercept = 100, linetype = 2,
                        linewidth = 0.4, color = "grey50") +
    ggplot2::geom_line(
      data = d0,
      ggplot2::aes(x = .data$time, y = .data$pct_smooth,
                   group = .data$cell_name),
      color = "grey75", linewidth = 0.4, alpha = 0.35, na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = d_mean,
      ggplot2::aes(x = .data$time, y = .data$mean_pct),
      color = "black", linewidth = 1.0, na.rm = TRUE
    ) +
    ggplot2::labs(x = "Time (s)", y = "% baseline") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
