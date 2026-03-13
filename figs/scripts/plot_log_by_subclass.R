library(dplyr)
library(ggplot2)

plot_log2fc_by_subclass <- function(df0,
                                    subclasses = c("L5_IT", "L5_ET"),
                                    x_bin_s = 1,                # bin width in seconds
                                    eps_pct = 0.5,              # floor for 0% baseline
                                    xlim = c(-10, 50),
                                    show_mean = TRUE,
                                    mean_lwd = 1.2,
                                    trace_lwd = 0.4,
                                    trace_alpha = 0.25,
                                    trace_color = "grey70",
                                    plot_title = NULL) {

  d <- df0 %>%
    dplyr::filter(assigned_subclass %in% subclasses) %>%
    dplyr::mutate(
      x_bin = floor(time / x_bin_s) * x_bin_s
    ) %>%
    dplyr::group_by(x_bin, cell_name, assigned_subclass, Species) %>%
    dplyr::summarise(
      y_pct = mean(percent_change, na.rm = TRUE),   # expects % baseline scale (100=baseline)
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      log2_fc = log2(pmax(y_pct, eps_pct) / 100)
    ) %>%
    dplyr::filter(x_bin >= xlim[1], x_bin <= xlim[2])

  # mean trace per facet (optionally by Species too)
  d_mean <- d %>%
    dplyr::group_by(assigned_subclass, x_bin) %>%
    dplyr::summarise(mean_log2_fc = mean(log2_fc, na.rm = TRUE), .groups = "drop")

  p <- ggplot(d, aes(x = x_bin, y = log2_fc, group = cell_name)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    facet_wrap(~assigned_subclass, nrow = 1) +
    coord_cartesian(xlim = xlim) +
    labs(
      title = plot_title,
      x = "Time (s)",
      y = expression(log[2]~"fold-change (baseline = 0)")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  if (isTRUE(show_mean)) {
    p <- p +
      geom_line(
        data = d_mean,
        aes(x = x_bin, y = mean_log2_fc, group = assigned_subclass),
        inherit.aes = FALSE,
        linewidth = mean_lwd
      )
  }

  p
}

# usage
p <- plot_log2fc_by_subclass(
  df0,
  subclasses = c("L5_ET"),
  x_bin_s = 1,
  eps_pct = 0.5,
  xlim = c(-10, 35),
  trace_alpha = 0.40,
  trace_color = "purple3",
  mean_lwd = 1.4
)


pKet <- df_ket %>%
  plot_log2fc_by_subclass(
  .,
  subclasses = c("L5_IT"),
  x_bin_s = 1,
  eps_pct = 0.5,
  xlim = c(-10, 35),
  trace_alpha = 0.40,
  trace_color = "purple3",
  mean_lwd = 1.4
)
print(pKet)

table(df_kcon$bath)
table(df_ket$bath)

df_kcon %>%
  #dplyr::filter(bath%in% "Ket") %>%
  dplyr::filter(bath == "none") %>%
  plot_log2fc_by_subclass(
  .,
  subclasses = c("L5_IT"),
  x_bin_s = 1,
  eps_pct = 0.5,
  xlim = c(-10, 35),
  trace_alpha = 0.40,
  trace_color = "purple3",
  mean_lwd = 1.4
)
print(pKet)

pKetIT <- df_ket %>%
  plot_log2fc_by_subclass(
  .,
  subclasses = c("L5_IT"),
  x_bin_s = 1,
  eps_pct = 0.5,
  xlim = c(-10, 35),
  trace_alpha = 0.40,
  trace_color = "purple3",
  mean_lwd = 1.4
)
print(pKet)



plot_log10fc_by_subclass <- function(df0,
                                     subclasses = c("L5_IT", "L5_ET"),
                                     x_bin_s = 1,
                                     eps_pct = 0.5,
                                     xlim = c(-10, 50),
                                     ylim = c(-1, 0.6),

                                     show_mean = TRUE,

                                     trace_lwd = 0.4,
                                     trace_alpha = 0.25,
                                     trace_color = "grey70",

                                     mean_lwd = 1.3,
                                     mean_color = "black",

                                     plot_title = NULL) {

  d <- df0 %>%
    dplyr::filter(assigned_subclass %in% subclasses) %>%
    dplyr::mutate(
      x_bin = floor(time / x_bin_s) * x_bin_s
    ) %>%
    dplyr::group_by(x_bin, cell_name, assigned_subclass, Species) %>%
    dplyr::summarise(
      y_pct = mean(percent_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cell_name, assigned_subclass) %>%
    tidyr::complete(
      x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)
    ) %>%
    dplyr::mutate(y_pct = dplyr::coalesce(y_pct, 0)) %>%
    dplyr::mutate(
      log10_fc = log10(pmax(y_pct, eps_pct) / 100)
    ) %>%
    dplyr::filter(x_bin >= xlim[1], x_bin <= xlim[2])

  d_mean <- d %>%
    dplyr::group_by(assigned_subclass, x_bin) %>%
    dplyr::summarise(mean_log10_fc = mean(log10_fc, na.rm = TRUE), .groups = "drop")

  ggplot(d, aes(x = x_bin, y = log10_fc, group = cell_name)) +

    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey40") +

    geom_line(color = trace_color,
              alpha = trace_alpha,
              linewidth = trace_lwd) +

    {
      if (show_mean)
        geom_line(data = d_mean,
                  aes(x = x_bin, y = mean_log10_fc),
                  inherit.aes = FALSE,
                  color = mean_color,
                  linewidth = mean_lwd)
    } +

    facet_wrap(~assigned_subclass, nrow = 1) +

    coord_cartesian(xlim = xlim, ylim = ylim) +

    scale_y_continuous(
      breaks = scales::pretty_breaks(n = 6)
    ) +

    labs(
      title = plot_title,
      x = "Time (s)",
      y = expression(log[10]~"fold-change from baseline")
    ) +

    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )
}

pET <- plot_log2fc_by_subclass(
  df0,
  subclasses = c("L5_ET"),
  x_bin_s = 1,
  trace_alpha = 0.4,
  trace_color = "purple1",
  xlim = c(-10,35)
  )
plot_log10fc_by_subclass(
  df_ket,
  subclasses = c("L5_ET"),
  x_bin_s = 1,
  trace_alpha = 0.4,
  trace_color = "purple1",
  ylim = c(-2, 0.8),
  xlim = c(-10,35),
  plot_title = "5HT response L5_ETs in Ketanserin"
)
pIT5 <- plot_log10fc_by_subclass(
  df0,
  subclasses = c("L5_IT"),
  x_bin_s = 1,
  trace_alpha = 0.4,
  trace_color = "purple1",
  ylim = c(-1, 0.8),
  xlim = c(-10,35),
  plot_title = "5HT response in L5_ITs"
)

print(p)
p23 <- plot_log10fc_by_subclass(
  df0,
  subclasses = c("L23_IT"),
  x_bin_s = 2,
  trace_alpha = 0.4,
  trace_color = "purple1",
  mean_color = "navyblue",
  ylim = c(-1, 0.8),
  xlim = c(-10,35)
)
p23 + ggtitle("5HT responses in L2/3 IT neurons")

plot_log10fc_by_subclass_depth <- function(df0,
                                           subclass_use = "L23_IT",
                                           depth_break = 500,
                                           x_bin_s = 1,
                                           eps_pct = 0.5,
                                           xlim = c(-10, 50),
                                           ylim = c(-1, 0.6),
                                           y_lims_by_row = NULL,
                                           row_order = c(1, 2),  # NEW: index-based row order (c(2,1) swaps)

                                           show_mean = TRUE,

                                           trace_lwd = 0.4,
                                           trace_alpha = 0.25,
                                           trace_color = "purple2",

                                           mean_lwd = 1.3,
                                           mean_color = "navyblue") {

  low_lbl  <- paste0("0-", depth_break)
  high_lbl <- paste0(">", depth_break)

  d <- df0 %>%
    dplyr::filter(assigned_subclass %in% subclass_use) %>%
    dplyr::mutate(
      depth_group = dplyr::case_when(
        !is.na(assigned_depth) & assigned_depth <= depth_break ~ low_lbl,
        !is.na(assigned_depth) & assigned_depth >  depth_break ~ high_lbl,
        TRUE ~ NA_character_
      ),
      x_bin = floor(time / x_bin_s) * x_bin_s
    ) %>%
    dplyr::filter(!is.na(depth_group)) %>%
    dplyr::group_by(x_bin, cell_name, assigned_subclass, Species, depth_group) %>%
    dplyr::summarise(
      y_pct = mean(percent_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cell_name, assigned_subclass, Species, depth_group) %>%
    tidyr::complete(
      x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)
    ) %>%
    dplyr::mutate(y_pct = dplyr::coalesce(y_pct, 0)) %>%
    dplyr::mutate(
      log10_fc = log10(pmax(y_pct, eps_pct) / 100)
    ) %>%
    dplyr::filter(x_bin >= xlim[1], x_bin <= xlim[2])

  # --- NEW: enforce stable + swappable facet row ordering (label-agnostic) ---
  stopifnot(is.numeric(row_order), length(row_order) == 2)
  stopifnot(setequal(as.integer(row_order), c(1L, 2L)))

  d <- d %>%
    dplyr::mutate(depth_group = factor(depth_group, levels = c(low_lbl, high_lbl)))

  depth_levels <- levels(d$depth_group)[as.integer(row_order)]
  d <- d %>%
    dplyr::mutate(depth_group = factor(depth_group, levels = depth_levels))

  row_levels <- levels(d$depth_group)
  # --------------------------------------------------------------------------

  d_mean <- d %>%
    dplyr::group_by(assigned_subclass, depth_group, x_bin) %>%
    dplyr::summarise(mean_log10_fc = mean(log10_fc, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(d, ggplot2::aes(x = x_bin, y = log10_fc, group = cell_name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    {
      if (show_mean)
        ggplot2::geom_line(
          data = d_mean,
          ggplot2::aes(x = x_bin, y = mean_log10_fc, group = 1),
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
    ggplot2::labs(
      title = paste0("By depth"),
      x = "Time (s)",
      y = expression(log[10]~"fold-change from baseline")
    ) +
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
          ggplot2::scale_y_continuous(
            limits = !!lim,
            breaks = scales::pretty_breaks(n = 6)
          )
        )
      )
    })

    p <- p + ggh4x::facetted_pos_scales(y = y_specs)
  }

  p
}


pDepth = plot_log10fc_by_subclass_depth(df0, trace_alpha = 0.4, depth_break = 400, x_bin_s = 2, y_lims_by_row = list(c(-.5, .6), c(-2, 0.55)), row_order = c(2,1))





plot_log10fc_by_subclass_cluster <- function(df0,
                                             subclass_use = "L23_IT",
                                             cluster_col = "cluster_Corr",
                                             cluster_levels = NULL,
                                             clusters_keep = NULL,
                                             x_bin_s = 1,
                                             eps_pct = 0.5,
                                             xlim = c(-10, 50),
                                             ylim = c(-1, 0.6),
                                             y_lims_by_row = NULL,
                                             show_mean = TRUE,
                                             trace_lwd = 0.4,
                                             trace_alpha = 0.25,
                                             trace_color = "purple2",
                                             mean_lwd = 1.3,
                                             mean_color = "navyblue") {

  stopifnot(cluster_col %in% names(df0))

  d <- df0 %>%
    dplyr::filter(.data$assigned_subclass %in% subclass_use) %>%
    dplyr::mutate(
      # IMPORTANT: keep as factor for facet ordering
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
    dplyr::summarise(
      y_pct = mean(.data$percent_change, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(cell_name, assigned_subclass, Species, cluster_group) %>%
    tidyr::complete(
      x_bin = seq(from = xlim[1], to = xlim[2], by = x_bin_s)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # interpret missing bins as "no spiking / no response"
      y_pct = dplyr::coalesce(y_pct, 0)
    ) %>%
    dplyr::mutate(
      log10_fc = log10(pmax(.data$y_pct, eps_pct) / 100)
    ) %>%
    dplyr::filter(.data$x_bin >= xlim[1], .data$x_bin <= xlim[2])

  d_mean <- d %>%
    dplyr::group_by(.data$assigned_subclass, .data$cluster_group, .data$x_bin) %>%
    dplyr::summarise(mean_log10_fc = mean(.data$log10_fc, na.rm = TRUE), .groups = "drop")

  # ... after d and d_mean are created ...

  # Determine the facet row levels actually present (in order)
  row_levels <- levels(d$cluster_group)
  row_levels <- row_levels[row_levels %in% unique(as.character(d$cluster_group))]

  # Build plot
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x_bin, y = .data$log10_fc, group = .data$cell_name)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, linewidth = 0.4, color = "grey40") +
    ggplot2::geom_line(color = trace_color, alpha = trace_alpha, linewidth = trace_lwd) +
    {
      if (show_mean)
        ggplot2::geom_line(
          data = d_mean,
          ggplot2::aes(x = .data$x_bin, y = .data$mean_log10_fc, group = 1),
          inherit.aes = FALSE,
          color = mean_color,
          linewidth = mean_lwd
        )
    } +
    ggplot2::facet_grid(cluster_group ~ assigned_subclass, scales = "free_y", drop = FALSE) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
    ggplot2::labs(
      #title = paste0("5-HT response in ", paste(subclass_use, collapse = ", ")," split by ", cluster_col),
      title = "Split by Cluster",
      x = "Time (s)",
      y = expression(log[10]~"fold-change from baseline")
    ) +
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

    # Create a list of formulas like: cluster_group == "<level>" ~ scale_y_continuous(limits=...)
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




# keep only two clusters
pCluster = plot_log10fc_by_subclass_cluster(
  df0,
  subclass_use = "L23_IT",
  cluster_levels = c("L2/3_IT_3", "L2/3_IT_1"),
  y_lims_by_row = list(c(-0.5, .6), c(-2, 0.55)),
  trace_alpha = 0.4,
  x_bin_s = 2,
  xlim = c(-10,35)
)


pDepth + ggtitle("Depth") | pCluster + ggtitle("Cluster")

