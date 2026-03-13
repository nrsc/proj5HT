plot_response_types_4panel_balanced_grey_spline <- function(df_prepped,
                                                            per_cell,
                                                            xlim = c(-5, 50),
                                                            k_smooth = 3,
                                                            # make traces faint
                                                            trace_alpha = 0.04,
                                                            trace_lwd = 0.25,
                                                            # baseline-centered transform
                                                            center = 100,
                                                            sigma = 6,
                                                            breaks = c(50, 70, 80, 90, 100, 110, 120, 140, 200),
                                                            # summary choices
                                                            summary_fun = c("median","mean"),
                                                            spline_df = 12,
                                                            spar = NULL) {

  summary_fun <- match.arg(summary_fun)

  keep_tbl <- per_cell %>%
    dplyr::mutate(
      resp_class = dplyr::case_when(
        response_type == "excitation" ~ "Excitation",
        response_type == "inhibition" ~ "Inhibition",
        response_type == "biphasic" & biphasic_order == "exc_then_inh" ~ "Biphasic exc→inh",
        response_type == "biphasic" & biphasic_order == "inh_then_exc" ~ "Biphasic inh→exc",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(resp_class)) %>%
    dplyr::select(cell_name, resp_class) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      resp_class = factor(resp_class,
                          levels = c("Excitation","Inhibition","Biphasic exc→inh","Biphasic inh→exc"))
    )

  # per-cell smoothed traces
  d <- df_prepped %>%
    dplyr::inner_join(keep_tbl, by = "cell_name") %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(resp_class, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth,
                                 fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  dt <- 0.2  # 0.1–1.0 depending on sampling; start with 0.2s
  d2 <- d %>%
    dplyr::mutate(tbin = dt * round(time / dt))

  # robust per-time summary (median usually better for skewed excitation)
  summariser <- if (summary_fun == "median") stats::median else base::mean

  summarize_on_grid <- function(d, dt = 0.1, fun = c("median","mean")) {
    fun <- match.arg(fun)
    f <- if (fun == "median") stats::median else base::mean

    grid <- seq(min(d$time, na.rm = TRUE), max(d$time, na.rm = TRUE), by = dt)

    d %>%
      dplyr::group_by(resp_class, cell_name) %>%
      dplyr::group_modify(~{
        dd <- .x %>% dplyr::filter(is.finite(time), is.finite(pct_smooth))
        if (nrow(dd) < 2) return(tibble::tibble(time = grid, y = NA_real_))
        # linear interpolation onto grid
        y <- stats::approx(dd$time, dd$pct_smooth, xout = grid, rule = 1)$y
        tibble::tibble(time = grid, y = y)
      }) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(resp_class, time) %>%
      dplyr::summarise(y = f(y, na.rm = TRUE), .groups = "drop")
  }

  d_sum <- summarize_on_grid(d, dt = 0.1, fun = summary_fun)


  # smoothing spline for the black line per panel
  # d_spline <- d_sum %>%
  #   dplyr::group_by(resp_class) %>%
  #   dplyr::group_modify(~{
  #     dd <- .x
  #     if (nrow(dd) < 6) return(dd %>% dplyr::mutate(y_spline = y))
  #     fit <- stats::smooth.spline(x = dd$time, y = dd$y, df = spline_df)
  #     dd$y_spline <- stats::predict(fit, x = dd$time)$y
  #     dd
  #   }) %>%
  #   dplyr::ungroup()

  min_cells_per_time <- 5  # now this will actually work

  d_sum <- d2 %>%
    dplyr::group_by(resp_class, tbin) %>%
    dplyr::summarise(
      n = sum(is.finite(pct_smooth)),
      y = stats::median(pct_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      y = dplyr::if_else(n >= min_cells_per_time, y, NA_real_)
    ) %>%
    dplyr::rename(time = tbin)


  # Fit splines only on timepoints with enough data
  d_spline <- d_sum %>%
    dplyr::group_by(resp_class) %>%
    dplyr::group_modify(~{
      dd <- .x %>% dplyr::filter(is.finite(time), is.finite(y))
      if (nrow(dd) < 6) return(.x %>% dplyr::mutate(y_spline = NA_real_))
      fit <- stats::smooth.spline(x = dd$time, y = dd$y, df = spline_df)
      dd$y_spline <- stats::predict(fit, x = dd$time)$y
      dd
    }) %>%
    dplyr::ungroup()




  ggplot(d, aes(x = time, y = pct_smooth, group = cell_name)) +
    facet_wrap(~resp_class, ncol = 2) +
    geom_hline(yintercept = center, linetype = 2, linewidth = 0.5, color = "grey40") +

    # very faint grey traces
    geom_line(color = "grey60", alpha = trace_alpha, linewidth = trace_lwd, na.rm = TRUE) +

    # black spline-smoothed summary curve
    geom_line(
      data = d_spline,
      aes(x = time, y = y_spline, group = resp_class),
      color = "black",
      linewidth = 1.6,
      inherit.aes = FALSE
    ) +

    scale_y_continuous(
      trans = centered_pseudolog_trans(center = center, sigma = sigma),
      breaks = breaks,
      labels = function(x) sprintf("%.0f", x)
    ) +
    labs(
      title = "Response types",
      # subtitle = paste0("Grey = per-cell smoothed traces (k=", k_smooth, "); black = ",
      #                   summary_fun, " across cells then smoothing spline (df=", spline_df, ")."),
      x = "Time (s)",
      y = "% baseline (log-balanced scale)"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

p <- plot_response_types_4panel_balanced_grey_spline(
  df_prepped, per_cell2,
  xlim = c(-5, 50),
  k_smooth = 3,
  trace_alpha = 0.3,   # start very low
  trace_lwd = 0.5,
  sigma = 6,
  summary_fun = "mean",
  spline_df = 48,
  spar = 0.3
)
print(p)

