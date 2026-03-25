# puff_plotting.R

#source("puff_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tibble)
})

response_cols <- c(
  strong_inhibition   = "#B2182B",
  moderate_inhibition = "#D6604D",
  mild_inhibition     = "#F4A582",
  no_change           = "grey80",
  mild_excitation     = "#92C5DE",
  moderate_excitation = "#4393C3",
  strong_excitation   = "#2166AC"
)

pretty_labels <- c(
  strong_inhibition   = "Strong inhibition",
  moderate_inhibition = "Moderate inhibition",
  mild_inhibition     = "Mild inhibition",
  no_change           = "No change",
  mild_excitation     = "Mild excitation",
  moderate_excitation = "Moderate excitation",
  strong_excitation   = "Strong excitation"
)

plot_cell_qc <- function(df_prepped,
                         cell_id,
                         smooth_method = "rollmean",
                         smooth_k = 3,
                         loess_span = 0.2,
                         rapid_window = c(0, 10),
                         mid_window = c(10, 40),
                         late_window = c(40, 60),
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
    ggplot2::annotate("rect", xmin = rapid_window[1], xmax = rapid_window[2],
                      ymin = -Inf, ymax = Inf, alpha = 0.06) +
    ggplot2::annotate("rect", xmin = mid_window[1], xmax = mid_window[2],
                      ymin = -Inf, ymax = Inf, alpha = 0.04) +
    ggplot2::annotate("rect", xmin = late_window[1], xmax = late_window[2],
                      ymin = -Inf, ymax = Inf, alpha = 0.05) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_hline(yintercept = 100 + min_excursion_pct, linetype = 3, linewidth = 0.35, color = "#2166AC") +
    ggplot2::geom_hline(yintercept = 100 - min_excursion_pct, linetype = 3, linewidth = 0.35, color = "#B2182B") +
    ggplot2::geom_line(ggplot2::aes(y = .data$pct_of_baseline), color = "grey70", linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = .data$pct_smooth), color = "black", linewidth = 0.75) +
    ggplot2::labs(
      title = cell_id,
      subtitle = first_non_missing(d$assigned_subclass),
      x = "Time (s)",
      y = "% baseline"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}

plot_response_composition <- function(per_cell,
                                      x = "assigned_subclass") {

  assert_has_cols(per_cell, c(x, "response_call"))

  per_cell %>%
    dplyr::filter(!is.na(.data$response_call)) %>%
    dplyr::count(.data[[x]], .data$response_call, name = "n") %>%
    dplyr::group_by(.data[[x]]) %>%
    dplyr::mutate(frac = .data$n / sum(.data$n)) %>%
    dplyr::ungroup() %>%
    ggplot2::ggplot(ggplot2::aes(x = .data[[x]], y = .data$frac, fill = .data$response_call)) +
    ggplot2::geom_col(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = response_cols, labels = pretty_labels, drop = FALSE) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      x = NULL,
      y = "Fraction of cells",
      fill = "Response class"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

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
    ggplot2::ggplot(ggplot2::aes(x = .data[[x]], y = .data[[y]], fill = .data$frac)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(.data$frac, accuracy = 1)), size = 3) +
    ggplot2::scale_fill_gradient(low = "white", high = "steelblue") +
    ggplot2::labs(x = NULL, y = NULL, fill = "Fraction") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
}

plot_centered_violin <- function(per_cell) {
  assert_has_cols(per_cell, c("assigned_subclass", "response_peak", "direction"))

  ggplot2::ggplot(
    per_cell %>% dplyr::filter(is.finite(.data$response_peak)),
    ggplot2::aes(x = .data$assigned_subclass, y = .data$response_peak)
  ) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4) +
    ggplot2::geom_violin(fill = "grey85", color = "black", linewidth = 0.35, trim = TRUE) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$direction),
      position = ggplot2::position_jitter(width = 0.12),
      size = 1.2,
      alpha = 0.55
    ) +
    ggplot2::scale_color_manual(
      values = c(
        inhibition = "#B2182B",
        excitation = "#2166AC",
        biphasic = "purple",
        no_change = "grey50"
      ),
      guide = "none"
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      x = NULL,
      y = "Peak % baseline"
    ) +
    ggplot2::theme_minimal(base_size = 12)
}

plot_exc_inh_traces <- function(df_prepped,
                                per_cell,
                                subclasses = c("L23_IT", "L5_ET", "L5_IT"),
                                xlim = c(-5, 60),
                                smooth_k = 3) {

  assert_has_cols(df_prepped, c("cell_name", "time", "pct_of_baseline"))
  assert_has_cols(per_cell, c("cell_name", "assigned_subclass", "response_type"))

  d0 <- df_prepped %>%
    # drop subclass from df_prepped so the join does not create .x/.y duplicates
    dplyr::select(-dplyr::any_of(c("assigned_subclass", "response_type"))) %>%
    dplyr::inner_join(
      per_cell %>%
        dplyr::select("cell_name", "assigned_subclass", "response_type"),
      by = "cell_name"
    ) %>%
    dplyr::filter(
      .data$assigned_subclass %in% subclasses,
      .data$response_type %in% c("excitation", "inhibition"),
      .data$time >= xlim[1],
      .data$time <= xlim[2]
    ) %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(.data$pct_of_baseline, k = smooth_k, fill = NA, align = "center")
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
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey50") +
    ggplot2::geom_line(
      data = d0,
      ggplot2::aes(x = .data$time, y = .data$pct_smooth, group = .data$cell_name),
      color = "grey75",
      linewidth = 0.4,
      alpha = 0.35,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = d_mean,
      ggplot2::aes(x = .data$time, y = .data$mean_pct),
      color = "black",
      linewidth = 1.0,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "% baseline"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
