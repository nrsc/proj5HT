library(dplyr)
library(ggplot2)
library(zoo)

make_plot_df <- function(df_prepped, per_cell2, xlim = c(-5, 56), k_smooth = 3) {

  keep_tbl <- per_cell2 %>%
    dplyr::select(cell_name, response_type, assigned_subclass) %>%
    dplyr::distinct()

  d <- df_prepped %>%
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%  # avoid .x/.y collision
    dplyr::inner_join(keep_tbl, by = "cell_name") %>%
    dplyr::mutate(
      time = as.numeric(time),
      pct_of_baseline = as.numeric(pct_of_baseline),
      response_type = as.character(response_type),
      assigned_subclass = as.character(assigned_subclass)
    ) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(response_type, assigned_subclass, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(response_type), !is.na(assigned_subclass))

  if (nrow(d) == 0) {
    stop("make_plot_df(): joined+filtered plot df is empty. Check xlim and that cell_name overlaps between df_prepped and per_cell2.")
  }

  d
}

balanced_log_trans <- function(center = 100, base = 10, scale = 1) {
  scales::trans_new(
    name = paste0("balanced_log_center_", center),
    transform = function(x) {
      d <- x - center
      sign(d) * (log(abs(d) / scale + 1, base = base))
    },
    inverse = function(y) {
      center + sign(y) * (base^(abs(y)) - 1) * scale
    }
  )
}


plot_traces_faceted_by_type_and_subclass <- function(
    df_prepped,
    per_cell,
    xlim = c(-5, 56),
    k_smooth = 3,
    alpha_traces = 0.20,
    trace_color = "grey70",
    mean_color = "black",
    mean_lwd = 1.1,

    # mean smoothing control
    mean_smooth = c("spline", "rollmean", "none"),
    spline_df = 20,
    mean_roll_k = 11,

    # balanced log axis
    y_bal_log = TRUE,
    y_center = 100,
    y_base = 10,
    y_scale = 5,   # bigger -> less compression near 100
    y_breaks = c(60, 70, 80, 90, 100, 110, 120, 140, 170, 200, 250)
) {
  mean_smooth <- match.arg(mean_smooth)

  # ---- keep only plotted classes (do NOT modify per_cell itself) ----
  keep_tbl <- per_cell %>%
    dplyr::mutate(response_type = as.character(response_type)) %>%
    dplyr::filter(!is.na(response_type), response_type != "no_change") %>%
    dplyr::select(cell_name, response_type, assigned_subclass) %>%
    dplyr::distinct()

  # ---- join + window ----
  d <- df_prepped %>%
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%
    dplyr::inner_join(keep_tbl, by = "cell_name") %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(response_type, assigned_subclass, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  if (nrow(d) == 0) stop("No rows to plot after join/xlim filter.")

  # ---- mean per time (separate data frame; avoids recycle/case_when issues) ----
  mean_df <- d %>%
    dplyr::filter(is.finite(pct_smooth), is.finite(time)) %>%
    dplyr::group_by(response_type, assigned_subclass, time) %>%
    dplyr::summarise(
      n_cells = dplyr::n_distinct(cell_name),
      mean_y = mean(pct_smooth, na.rm = TRUE),
      .groups = "drop"
    )

  # optional smoothing of mean curve
  mean_df <- mean_df %>%
    dplyr::group_by(response_type, assigned_subclass) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      mean_y_smooth =
        dplyr::case_when(
          mean_smooth == "none" ~ mean_y,
          mean_smooth == "rollmean" ~ as.numeric(zoo::rollmean(mean_y, k = mean_roll_k, fill = NA, align = "center")),
          mean_smooth == "spline" ~ {
            ok <- is.finite(mean_y) & is.finite(time)
            if (sum(ok) < 6) rep(NA_real_, length(mean_y)) else {
              fit <- stats::smooth.spline(x = time[ok], y = mean_y[ok], df = spline_df)
              out <- rep(NA_real_, length(mean_y))
              out[ok] <- stats::predict(fit, x = time[ok])$y
              out
            }
          },
          TRUE ~ mean_y
        )
    ) %>%
    dplyr::ungroup()

  # ---- plot ----
  p <- ggplot2::ggplot(d, ggplot2::aes(x = time, y = pct_smooth, group = cell_name)) +
    ggplot2::facet_grid(response_type ~ assigned_subclass, drop = FALSE) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey50") +
    ggplot2::geom_line(color = trace_color, alpha = alpha_traces, linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(
      data = mean_df,
      ggplot2::aes(x = time, y = mean_y_smooth),
      inherit.aes = FALSE,
      color = mean_color,
      linewidth = mean_lwd,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "% baseline",
      title = "Puff response traces by response type × subclass",
      subtitle = paste0("Traces smoothed (k=", k_smooth, "); mean_smooth=", mean_smooth,
                        if (mean_smooth == "spline") paste0(" (df=", spline_df, ")") else "")
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(y_bal_log)) {
    p <- p + ggplot2::scale_y_continuous(
      trans = balanced_log_trans(center = y_center, base = y_base, scale = y_scale),
      breaks = y_breaks,
      labels = y_breaks
    )
  }

  p
}

p <- plot_traces_faceted_by_type_and_subclass(
  df_prepped, per_cell2,
  xlim = c(-5, 50),
  mean_smooth = "spline",
  #mean_roll_k = 3,
  k_smooth = 3,
  spline_df = 48,
  trace_color = "grey75",
  alpha_traces = 0.5,
  mean_color = "black",
  y_scale = 20
)
print(p)

p <- plot_traces_faceted_by_type_and_subclass(
  df_prepped, per_cell2,
  mean_smooth = "none",
  trace_color = "grey75",
  mean_color = "white",
  alpha_traces = 0.2
)
print(p)



p_mean <- plot_traces_faceted_by_type_and_subclass(df_prepped, per_cell2,
                                                   xlim = c(-5, 56),
                                                   k_smooth = 3,
                                                   min_cells_per_time = 1,
                                                   mean_smooth = "spline",
                                                   spline_df = 1)
print(p_mean)


