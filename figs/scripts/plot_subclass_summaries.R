# =========================================================
# Plot functions (updated)
# - Explicit argument names (per_cell / per_cell2)
# - Safer joins (no accidental column duplication)
# - Optional facet scaling control for debugging
# =========================================================

balanced_log_trans <- function(center = 100, base = 10, scale = 5) {
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

balanced_log100_trans <- function(scale = 10) {
  scales::trans_new(
    name = paste0("blog100_", scale),
    transform = function(y) {
      d <- y - 100
      100 + sign(d) * log1p(abs(d) / scale)
    },
    inverse = function(z) {
      d <- z - 100
      100 + sign(d) * scale * (exp(abs(d)) - 1)
    }
  )
}

make_plot_df <- function(df_prepped,
                         per_cell,
                         xlim = c(-5, 56),
                         k_smooth = 3,
                         drop_no_change = TRUE) {

  stopifnot(is.data.frame(df_prepped), is.data.frame(per_cell))
  need_df  <- c("cell_name", "time", "pct_of_baseline")
  need_pc  <- c("cell_name", "assigned_subclass", "response_type")
  if (!all(need_df %in% names(df_prepped))) stop("df_prepped missing: ", paste(setdiff(need_df, names(df_prepped)), collapse=", "))
  if (!all(need_pc %in% names(per_cell)))   stop("per_cell missing: ", paste(setdiff(need_pc, names(per_cell)), collapse=", "))

  per_cell0 <- per_cell %>%
    dplyr::mutate(
      response_type = as.character(.data$response_type),
      biphasic_order = if ("biphasic_order" %in% names(per_cell)) as.character(.data$biphasic_order) else NA_character_
    ) %>%
    dplyr::filter(!is.na(.data$response_type)) %>%
    { if (drop_no_change) dplyr::filter(., .data$response_type != "no_change") else . } %>%
    # drop unknown order ONLY for biphasic rows
    dplyr::filter(!( .data$response_type == "biphasic" & is.na(.data$biphasic_order))) %>%
    dplyr::select(.data$cell_name, .data$response_type, .data$assigned_subclass, .data$biphasic_order) %>%
    dplyr::distinct()

  d <- df_prepped %>%
    # avoid duplicating subclass column if df_prepped also has it
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%
    dplyr::inner_join(per_cell0, by = "cell_name") %>%
    dplyr::mutate(time = as.numeric(.data$time)) %>%
    dplyr::filter(.data$time >= xlim[1], .data$time <= xlim[2]) %>%
    dplyr::group_by(.data$response_type, .data$assigned_subclass, .data$cell_name) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
      #pct_smooth = dplyr::case_when(
      #  !is.finite(pct_smooth) & is.finite(pct_of_baseline) ~ 0,  # true silent bins → 0
      #  TRUE ~ pct_smooth
      #)
    ) %>%
    dplyr::ungroup()

  # d <- df_prepped %>%
  #   dplyr::select(-dplyr::any_of("assigned_subclass")) %>%
  #   dplyr::inner_join(keep_tbl, by = "cell_name") %>%
  #   dplyr::mutate(time = as.numeric(time)) %>%
  #   dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
  #   dplyr::group_by(response_type, assigned_subclass, cell_name) %>%
  #   dplyr::arrange(time, .by_group = TRUE) %>%
  #   dplyr::mutate(
  #     pct0 = tidyr::replace_na(.data$pct_of_baseline, 0),
  #     pct_smooth = zoo::rollmean(.data$pct0, k = k_smooth, fill = NA_real_, align = "center"),
  #     pct_smooth = dplyr::if_else(!is.finite(.data$pct_smooth), .data$pct0, .data$pct_smooth)
  #   ) %>%
  #   dplyr::select(-.data$pct0) %>%
  #   dplyr::ungroup()


  if (nrow(d) == 0) stop("make_plot_df: no rows after join/xlim filters.")
  d
}

make_plot_df_ei <- function(df_prepped, per_cell,
                            xlim = c(-5, 56),
                            k_smooth = 3,
                            time_bin_s = NULL,          # if NULL, auto from dt_med_overall
                            min_cells_per_time = 3) {

  stopifnot(is.data.frame(df_prepped), is.data.frame(per_cell))

  if (is.null(time_bin_s)) {
    dt_est <- df_prepped %>%
      dplyr::group_by(.data$cell_name) %>%
      dplyr::summarise(dt = stats::median(diff(sort(unique(.data$time))), na.rm = TRUE), .groups="drop") %>%
      dplyr::summarise(dt = stats::median(.data$dt, na.rm = TRUE), .groups="drop") %>%
      dplyr::pull(.data$dt)
    time_bin_s <- max(0.1, dt_est)  # floor at 100 ms
  }

  d0 <- df_prepped %>%
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%
    dplyr::inner_join(
      per_cell %>% dplyr::select(.data$cell_name, .data$response_type, .data$assigned_subclass),
      by = "cell_name"
    ) %>%
    dplyr::filter(
      .data$response_type %in% c("excitation", "inhibition"),
      .data$time >= xlim[1], .data$time <= xlim[2],
      is.finite(.data$pct_of_baseline)
    ) %>%
    dplyr::mutate(time = as.numeric(.data$time))

  d_traces <- d0 %>%
    dplyr::group_by(.data$response_type, .data$assigned_subclass, .data$cell_name) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
      # pct_smooth = dplyr::case_when(
      #   !is.finite(pct_smooth) & is.finite(pct_of_baseline) ~ 0,  # true silent bins → 0
      #   TRUE ~ pct_smooth
      # )
    ) %>%
    dplyr::ungroup()

  # d_traces <- d0 %>%
  #   dplyr::group_by(response_type, assigned_subclass, cell_name) %>%
  #   dplyr::arrange(time, .by_group = TRUE) %>%
  #   dplyr::mutate(
  #     pct0 = tidyr::replace_na(.data$pct_of_baseline, 0),
  #     pct_smooth = zoo::rollmean(.data$pct0, k = k_smooth, fill = NA_real_, align = "center"),
  #     pct_smooth = dplyr::if_else(!is.finite(.data$pct_smooth), .data$pct0, .data$pct_smooth)
  #   ) %>%
  #   dplyr::select(-.data$pct0) %>%
  #   dplyr::ungroup()


  d_mean <- d_traces %>%
    dplyr::mutate(time_bin = round(.data$time / time_bin_s) * time_bin_s) %>%
    dplyr::group_by(.data$response_type, .data$assigned_subclass, .data$time_bin) %>%
    dplyr::summarise(
      n_cells = dplyr::n_distinct(.data$cell_name[is.finite(.data$pct_smooth)]),
      mean_pct = mean(.data$pct_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(.data$mean_pct), .data$n_cells >= min_cells_per_time) %>%
    dplyr::rename(time = .data$time_bin)

  list(traces = d_traces, mean = d_mean, time_bin_s = time_bin_s)
}

make_mean_df <- function(d,
                         mean_smooth = c("spline","rollmean","none"),
                         spline_df = 20,
                         mean_roll_k = 11) {

  mean_smooth <- match.arg(mean_smooth)

  mean_df <- d %>%
    dplyr::filter(is.finite(.data$pct_smooth), is.finite(.data$time)) %>%
    dplyr::group_by(.data$panel, .data$assigned_subclass, .data$time) %>%
    dplyr::summarise(
      n_cells = dplyr::n_distinct(.data$cell_name),
      mean_y = mean(.data$pct_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$panel, .data$assigned_subclass) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(
      mean_y_smooth = dplyr::case_when(
        mean_smooth == "none" ~ .data$mean_y,
        mean_smooth == "rollmean" ~ as.numeric(zoo::rollmean(.data$mean_y, k = mean_roll_k, fill = NA, align = "center")),
        mean_smooth == "spline" ~ {
          ok <- is.finite(.data$mean_y) & is.finite(.data$time)
          if (sum(ok) < 6) rep(NA_real_, length(.data$mean_y)) else {
            fit <- stats::smooth.spline(x = .data$time[ok], y = .data$mean_y[ok], df = spline_df)
            out <- rep(NA_real_, length(.data$mean_y))
            out[ok] <- stats::predict(fit, x = .data$time[ok])$y
            out
          }
        },
        TRUE ~ .data$mean_y
      )
    ) %>%
    dplyr::ungroup()

  mean_df
}

# ---------------------------------------------------------
# 1) Excitation vs inhibition (non-biphasic)
# ---------------------------------------------------------
plot_exc_inh_by_subclass <- function(df_prepped, per_cell,
                                     xlim = c(-5, 56),
                                     k_smooth = 3,
                                     min_cells_per_time = 3,
                                     time_bin_s = NULL,
                                     mean_smooth = c("none","spline"),
                                     spline_df = 25,
                                     trace_color = "grey80",
                                     mean_color = "black",
                                     alpha_traces = 0.35,
                                     lw_traces = 0.6,
                                     lw_mean = 1.2,
                                     y_scale = NULL,          # if not NULL => balanced log about 100
                                     facet_scales = c("free_y","fixed")) {

  mean_smooth <- match.arg(mean_smooth)
  facet_scales <- match.arg(facet_scales)

  dd <- make_plot_df_ei(df_prepped, per_cell,
                        xlim = xlim,
                        k_smooth = k_smooth,
                        time_bin_s = time_bin_s,
                        min_cells_per_time = min_cells_per_time)

  d_traces <- dd$traces
  d_mean   <- dd$mean

  if (mean_smooth == "spline" && nrow(d_mean) > 0) {
    d_mean <- d_mean %>%
      dplyr::group_by(.data$response_type, .data$assigned_subclass) %>%
      dplyr::arrange(.data$time, .by_group = TRUE) %>%
      dplyr::mutate(
        mean_pct = stats::smooth.spline(.data$time, .data$mean_pct, df = spline_df)$y
      ) %>%
      dplyr::ungroup()
  }

  p <- ggplot2::ggplot() +
    ggplot2::facet_grid(response_type ~ assigned_subclass, scales = facet_scales) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5, color = "grey50") +
    ggplot2::geom_line(
      data = d_traces,
      ggplot2::aes(x = time, y = pct_smooth, group = cell_name),
      color = trace_color,
      alpha = alpha_traces,
      linewidth = lw_traces,
      na.rm = TRUE
    ) +
    ggplot2::geom_line(
      data = d_mean,
      ggplot2::aes(x = time, y = mean_pct),
      color = mean_color,
      linewidth = lw_mean,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "% baseline",
      title = "Excitation vs inhibition by subclass",
      subtitle = "Non-biphasic responses"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(y_scale)) {
    p <- p + ggplot2::scale_y_continuous(trans = balanced_log100_trans(scale = y_scale))
  }

  p
}

# ---------------------------------------------------------
# 2) Biphasic modalities (exc→inh vs inh→exc)
# ---------------------------------------------------------
plot_biphasic_modalities_by_subclass <- function(df_prepped,
                                                 per_cell,
                                                 xlim = c(-5, 50),
                                                 k_smooth = 3,
                                                 alpha_traces = 0.20,
                                                 trace_color = "grey75",
                                                 mean_color = "black",
                                                 mean_lwd = 1.1,
                                                 mean_smooth = c("spline","rollmean","none"),
                                                 spline_df = 20,
                                                 mean_roll_k = 11,
                                                 y_bal_log = TRUE,
                                                 y_center = 100,
                                                 y_base = 10,
                                                 y_scale = 8,
                                                 y_breaks = c(60,70,80,90,100,110,120,140,170,200,250),
                                                 facet_scales = c("free_y","fixed")) {

  mean_smooth <- match.arg(mean_smooth)
  facet_scales <- match.arg(facet_scales)

  d <- make_plot_df(df_prepped, per_cell, xlim = xlim, k_smooth = k_smooth, drop_no_change = TRUE) %>%
    dplyr::filter(.data$response_type == "biphasic") %>%
    dplyr::filter(.data$biphasic_order %in% c("exc_then_inh", "inh_then_exc")) %>%  # drop unknown
    dplyr::mutate(
      panel = dplyr::case_when(
        .data$biphasic_order == "exc_then_inh" ~ "Biphasic exc→inh",
        .data$biphasic_order == "inh_then_exc" ~ "Biphasic inh→exc"
      ),
      panel = factor(.data$panel, levels = c("Biphasic exc→inh","Biphasic inh→exc"))
    )

  if (nrow(d) == 0) stop("plot_biphasic_modalities_by_subclass: no biphasic rows after filtering.")

  mean_df <- make_mean_df(
    d %>% dplyr::select(-dplyr::any_of("response_type")),
    mean_smooth = mean_smooth,
    spline_df = spline_df,
    mean_roll_k = mean_roll_k
  )

  p <- ggplot2::ggplot(d, ggplot2::aes(.data$time, .data$pct_smooth, group = .data$cell_name)) +
    ggplot2::facet_grid(panel ~ assigned_subclass, drop = FALSE, scales = facet_scales) +
    ggplot2::geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4, color = "grey50") +
    ggplot2::geom_line(color = trace_color, alpha = alpha_traces, linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_line(
      data = mean_df,
      ggplot2::aes(.data$time, .data$mean_y_smooth),
      inherit.aes = FALSE,
      color = mean_color,
      linewidth = mean_lwd,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "% baseline",
      title = "Biphasic modalities by subclass"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(y_bal_log)) {
    p <- p + ggplot2::scale_y_continuous(
      trans = balanced_log_trans(center = y_center, base = y_base, scale = y_scale),
      breaks = y_breaks, labels = y_breaks
    )
  }

  p
}


plot_4row_modalities_by_subclass <- function(
    df_prepped, per_cell,
    xlim = c(-5, 39),
    k_smooth = 11,
    min_cells_per_time = 2,
    time_bin_s = NULL,

    mean_smooth = c("spline","none"),
    spline_df = 30,

    # ---- TRACE APPEARANCE ----
    trace_alpha = 0.30,
    trace_color = "grey75",

    # optional biphasic override
    trace_alpha_bi = NULL,
    trace_color_bi = NULL,

    # ---- MEAN APPEARANCE ----
    mean_alpha = 1,
    mean_color = "black",
    mean_lwd = 1.2,

    facet_scales = c("fixed","free_y"),

    y_bal_log = FALSE,
    y_center = 100,
    y_base = 10,
    y_scale = 8,
    y_breaks = c(60,70,80,90,100,110,120,140,170,200,250)
) {

  mean_smooth <- match.arg(mean_smooth)
  facet_scales <- match.arg(facet_scales)

  # fallback so biphasic uses same styling unless specified
  if (is.null(trace_alpha_bi)) trace_alpha_bi <- trace_alpha
  if (is.null(trace_color_bi)) trace_color_bi <- trace_color

  # ---------------- EI ----------------
  dd_ei <- make_plot_df_ei(df_prepped, per_cell,
                           xlim = xlim,
                           k_smooth = k_smooth,
                           time_bin_s = time_bin_s,
                           min_cells_per_time = min_cells_per_time)

  d_ei_tr <- dd_ei$traces %>%
    dplyr::mutate(
      panel = dplyr::case_when(
        response_type == "excitation" ~ "Excitation",
        response_type == "inhibition" ~ "Inhibition"
      ),
      trace_alpha = trace_alpha,
      trace_color = trace_color
    ) %>%
    dplyr::select(panel, assigned_subclass, cell_name, time,
                  y = pct_smooth, trace_alpha, trace_color)

  d_ei_mean <- dd_ei$mean %>%
    dplyr::mutate(
      panel = dplyr::case_when(
        response_type == "excitation" ~ "Excitation",
        response_type == "inhibition" ~ "Inhibition"
      )
    ) %>%
    dplyr::select(panel, assigned_subclass, time, mean_y = mean_pct)

  if (mean_smooth == "spline" && nrow(d_ei_mean) > 0) {
    d_ei_mean <- d_ei_mean %>%
      dplyr::group_by(panel, assigned_subclass) %>%
      dplyr::arrange(time, .by_group = TRUE) %>%
      dplyr::mutate(mean_y = stats::smooth.spline(time, mean_y, df = spline_df)$y) %>%
      dplyr::ungroup()
  }

  # ---------------- BIPHASIC ----------------
  d_bi_base <- make_plot_df(df_prepped, per_cell,
                            xlim = xlim, k_smooth = k_smooth, drop_no_change = TRUE) %>%
    dplyr::filter(response_type == "biphasic",
                  biphasic_order %in% c("inh_then_exc","exc_then_inh"))

  d_bi_tr <- d_bi_base %>%
    dplyr::mutate(
      panel = dplyr::case_when(
        biphasic_order == "inh_then_exc" ~ "Biphasic inh→exc",
        biphasic_order == "exc_then_inh" ~ "Biphasic exc→inh"
      ),
      trace_alpha = trace_alpha_bi,
      trace_color = trace_color_bi
    ) %>%
    dplyr::select(panel, assigned_subclass, cell_name, time,
                  y = pct_smooth, trace_alpha, trace_color)

  mean_df_bi <- make_mean_df(
    d_bi_base %>% dplyr::select(panel = biphasic_order,
                                assigned_subclass, cell_name, time, pct_smooth),
    mean_smooth = if (mean_smooth == "spline") "spline" else "none",
    spline_df = spline_df
  ) %>%
    dplyr::mutate(
      panel = dplyr::recode(panel,
                            "inh_then_exc" = "Biphasic inh→exc",
                            "exc_then_inh" = "Biphasic exc→inh")
    ) %>%
    dplyr::select(panel, assigned_subclass, time, mean_y = mean_y_smooth)

  # ---------------- COMBINE ----------------
  d_tr   <- dplyr::bind_rows(d_ei_tr, d_bi_tr)
  d_mean <- dplyr::bind_rows(d_ei_mean, mean_df_bi)

  panel_levels <- c("Excitation","Inhibition",
                    "Biphasic inh→exc","Biphasic exc→inh")

  d_tr$panel   <- factor(d_tr$panel,   levels = panel_levels)
  d_mean$panel <- factor(d_mean$panel, levels = panel_levels)

  # ---------------- PLOT ----------------
  p <- ggplot2::ggplot() +
    ggplot2::facet_grid(panel ~ assigned_subclass,
                        scales = facet_scales, drop = FALSE) +

    ggplot2::geom_hline(yintercept = 100,
                        linetype = 2, linewidth = 0.4, color = "grey50") +

    ggplot2::geom_line(
      data = d_tr,
      ggplot2::aes(time, y, group = cell_name,
                   alpha = trace_alpha, colour = trace_color),
      linewidth = 0.5, na.rm = TRUE, show.legend = FALSE
    ) +

    ggplot2::geom_line(
      data = d_mean,
      ggplot2::aes(time, mean_y),
      colour = mean_color,
      alpha = mean_alpha,
      linewidth = mean_lwd,
      na.rm = TRUE
    ) +

    ggplot2::scale_alpha_identity() +
    ggplot2::scale_colour_identity() +

    ggplot2::labs(
      x = "Time (s)",
      y = "% baseline",
      title = "Response modalities by subclass (4-row)"
    ) +

    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (isTRUE(y_bal_log)) {
    p <- p + ggplot2::scale_y_continuous(
      trans = balanced_log_trans(center = y_center,
                                 base = y_base,
                                 scale = y_scale),
      breaks = y_breaks, labels = y_breaks
    )
  }

  p
}



