library(dplyr)
library(ggplot2)
library(tidyr)
library(scales)

plot_pct_by_subclass <- function(df0,
                                 subclasses = c("L5_IT", "L5_ET"),
                                 x_bin_s = 1,
                                 xlim = c(-10, 50),
                                 ylim = c(60, 140),   # around baseline=100
                                 show_mean = TRUE,
                                 trace_lwd = 0.4,
                                 trace_alpha = 0.25,
                                 trace_color = "grey70",
                                 mean_lwd = 1.3,
                                 mean_color = "black",
                                 plot_title = "5-HT response by subclass (% baseline)") {

  d <- df0 %>%
    dplyr::filter(.data$assigned_subclass %in% subclasses) %>%
    dplyr::mutate(x_bin = floor(.data$time / x_bin_s) * x_bin_s) %>%
    dplyr::group_by(.data$x_bin, .data$cell_name, .data$assigned_subclass, .data$Species) %>%
    dplyr::summarise(y_pct = mean(.data$percent_change, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$cell_name, .data$assigned_subclass, .data$Species) %>%
    tidyr::complete(x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$x_bin >= xlim[1], .data$x_bin <= xlim[2])

  d_mean <- d %>%
    dplyr::group_by(.data$assigned_subclass, .data$x_bin) %>%
    dplyr::summarise(mean_y_pct = mean(.data$y_pct, na.rm = TRUE), .groups = "drop")

  ggplot2::ggplot(d, ggplot2::aes(x = .data$x_bin, y = .data$y_pct, group = .data$cell_name)) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    {
      if (isTRUE(show_mean))
        ggplot2::geom_line(
          data = d_mean,
          ggplot2::aes(x = .data$x_bin, y = .data$mean_y_pct, group = 1),
          inherit.aes = FALSE,
          color = mean_color,
          linewidth = mean_lwd
        )
    } +
    ggplot2::facet_wrap(~assigned_subclass, nrow = 1) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::labs(title = plot_title, x = "Time (s)", y = "% of baseline (baseline = 100)") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )
}


plot_pct_by_subclass_depth <- function(df0,
                                       subclass_use = "L23_IT",
                                       depth_break = 500,
                                       x_bin_s = 1,
                                       xlim = c(-10, 50),
                                       ylim = c(60, 140),
                                       y_lims_by_row = NULL,
                                       row_order = c(1, 2),
                                       show_mean = TRUE,
                                       trace_lwd = 0.4,
                                       trace_alpha = 0.25,
                                       trace_color = "purple2",
                                       mean_lwd = 1.3,
                                       mean_color = "navyblue",
                                       plot_title = "By depth (% baseline)") {

  low_lbl  <- paste0("0-", depth_break)
  high_lbl <- paste0(">", depth_break)

  d <- df0 %>%
    dplyr::filter(.data$assigned_subclass %in% subclass_use) %>%
    dplyr::mutate(
      depth_group = dplyr::case_when(
        !is.na(.data$assigned_depth) & .data$assigned_depth <= depth_break ~ low_lbl,
        !is.na(.data$assigned_depth) & .data$assigned_depth >  depth_break ~ high_lbl,
        TRUE ~ NA_character_
      ),
      x_bin = floor(.data$time / x_bin_s) * x_bin_s
    ) %>%
    dplyr::filter(!is.na(.data$depth_group)) %>%
    dplyr::group_by(.data$x_bin, .data$cell_name, .data$assigned_subclass, .data$Species, .data$depth_group) %>%
    dplyr::summarise(y_pct = mean(.data$percent_change, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$cell_name, .data$assigned_subclass, .data$Species, .data$depth_group) %>%
    tidyr::complete(x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$x_bin >= xlim[1], .data$x_bin <= xlim[2])

  # stable + swappable ordering (same idea as your log10 fn)
  stopifnot(is.numeric(row_order), length(row_order) == 2)
  stopifnot(setequal(as.integer(row_order), c(1L, 2L)))

  d <- d %>%
    dplyr::mutate(depth_group = factor(.data$depth_group, levels = c(low_lbl, high_lbl)))

  depth_levels <- levels(d$depth_group)[as.integer(row_order)]
  d <- d %>%
    dplyr::mutate(depth_group = factor(.data$depth_group, levels = depth_levels))

  row_levels <- levels(d$depth_group)

  d_mean <- d %>%
    dplyr::group_by(.data$assigned_subclass, .data$depth_group, .data$x_bin) %>%
    dplyr::summarise(mean_y_pct = mean(.data$y_pct, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x_bin, y = .data$y_pct, group = .data$cell_name)) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    {
      if (show_mean)
        ggplot2::geom_line(
          data = d_mean,
          ggplot2::aes(x = .data$x_bin, y = .data$mean_y_pct, group = 1),
          inherit.aes = FALSE,
          color = mean_color,
          linewidth = mean_lwd
        )
    } +
    ggplot2::facet_grid(depth_group ~ assigned_subclass, scales = "free_y", drop = FALSE) +
    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = if (is.null(y_lims_by_row)) ylim else NULL
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::labs(title = plot_title, x = "Time (s)", y = "% of baseline (baseline = 100)") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(y_lims_by_row)) {
    stopifnot(length(y_lims_by_row) == length(row_levels))
    stopifnot(all(vapply(y_lims_by_row, function(z) is.numeric(z) && length(z) == 2, logical(1))))

    y_specs <- lapply(seq_along(row_levels), function(i) {
      lvl <- row_levels[[i]]
      lim <- y_lims_by_row[[i]]
      rlang::new_formula(
        lhs = rlang::expr(depth_group == !!lvl),
        rhs = rlang::expr(
          ggplot2::scale_y_continuous(limits = !!lim, breaks = scales::pretty_breaks(n = 6))
        )
      )
    })

    p <- p + ggh4x::facetted_pos_scales(y = y_specs)
  }

  p
}


plot_pct_by_subclass_cluster <- function(df0,
                                         subclass_use = "L23_IT",
                                         cluster_col = "cluster_Corr",
                                         cluster_levels = NULL,
                                         clusters_keep = NULL,
                                         x_bin_s = 1,
                                         xlim = c(-10, 50),
                                         ylim = c(60, 140),
                                         y_lims_by_row = NULL,
                                         show_mean = TRUE,
                                         trace_lwd = 0.4,
                                         trace_alpha = 0.25,
                                         trace_color = "purple2",
                                         mean_lwd = 1.3,
                                         mean_color = "navyblue",
                                         plot_title = "Split by Cluster (% baseline)") {

  stopifnot(cluster_col %in% names(df0))

  d <- df0 %>%
    dplyr::filter(.data$assigned_subclass %in% subclass_use) %>%
    dplyr::mutate(
      cluster_group = {
        x <- .data[[cluster_col]]
        if (!is.null(cluster_levels)) factor(x, levels = cluster_levels) else factor(x)
      },
      x_bin = floor(.data$time / x_bin_s) * x_bin_s
    ) %>%
    dplyr::filter(!is.na(.data$cluster_group)) %>%
    {
      if (!is.null(clusters_keep)) dplyr::filter(., as.character(.data$cluster_group) %in% clusters_keep) else .
    } %>%
    dplyr::group_by(.data$x_bin, .data$cell_name, .data$assigned_subclass, .data$Species, .data$cluster_group) %>%
    dplyr::summarise(y_pct = mean(.data$percent_change, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$cell_name, .data$assigned_subclass, .data$Species, .data$cluster_group) %>%
    tidyr::complete(x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)) %>%
    dplyr::mutate(
      # interpret missing bins as "no spiking / no response"
      y_pct = dplyr::coalesce(y_pct, 0)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$x_bin >= xlim[1], .data$x_bin <= xlim[2])

  d_mean <- d %>%
    dplyr::group_by(.data$assigned_subclass, .data$cluster_group, .data$x_bin) %>%
    dplyr::summarise(mean_y_pct = mean(.data$y_pct, na.rm = TRUE), .groups = "drop")

  row_levels <- levels(d$cluster_group)
  row_levels <- row_levels[row_levels %in% unique(as.character(d$cluster_group))]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x_bin, y = .data$y_pct, group = .data$cell_name)) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    {
      if (show_mean)
        ggplot2::geom_line(
          data = d_mean,
          ggplot2::aes(x = .data$x_bin, y = .data$mean_y_pct, group = 1),
          inherit.aes = FALSE,
          color = mean_color,
          linewidth = mean_lwd
        )
    } +
    ggplot2::facet_grid(cluster_group ~ assigned_subclass, scales = "free_y", drop = FALSE) +
    ggplot2::coord_cartesian(
      xlim = xlim,
      ylim = if (is.null(y_lims_by_row)) ylim else NULL
    ) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::labs(title = plot_title, x = "Time (s)", y = "% of baseline (baseline = 100)") +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold")
    )

  if (!is.null(y_lims_by_row)) {
    stopifnot(length(y_lims_by_row) == length(row_levels))
    stopifnot(all(vapply(y_lims_by_row, function(z) is.numeric(z) && length(z) == 2, logical(1))))

    y_specs <- lapply(seq_along(row_levels), function(i) {
      lvl <- row_levels[[i]]
      lim <- y_lims_by_row[[i]]
      rlang::new_formula(
        lhs = rlang::expr(cluster_group == !!lvl),
        rhs = rlang::expr(ggplot2::scale_y_continuous(limits = !!lim))
      )
    })

    p <- p + ggh4x::facetted_pos_scales(y = y_specs)
  }

  p
}


p_pct <- plot_pct_by_subclass(df0, subclasses = c("L5_IT","L5_ET"),
                              x_bin_s = 1, xlim = c(-10,35), ylim = c(0,300),
                              trace_alpha = 0.4, trace_color = "purple3")
p_pct

pDepth_pct <- plot_pct_by_subclass_depth(
  df0, subclass_use = "L23_IT",
  depth_break = 400, x_bin_s = 2, xlim = c(-10,35),
  y_lims_by_row = list(c(0, 300), c(0, 140)),
  row_order = c(2,1),
  trace_alpha = 0.4
)

pCluster_pct <- plot_pct_by_subclass_cluster(
  df0, subclass_use = "L23_IT",
  cluster_levels = c("L2/3_IT_1", "L2/3_IT_3"),
  x_bin_s = 2, xlim = c(-10,35),
  y_lims_by_row = list(c(80, 300), c(0, 300)),
  trace_alpha = 0.4
)

pDepth_pct + ggtitle("Depth (% baseline)") | pCluster_pct + ggtitle("Cluster (% baseline)")

