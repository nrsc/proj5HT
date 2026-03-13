suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(patchwork)
})

# ---------------------------------------------------------
# Helper: add smoothed trace (same spirit as your pipeline)
# ---------------------------------------------------------
add_smoothed_trace <- function(df_cell,
                               k_smooth = 11,
                               smooth_method = c("rollmean","rollmedian","loess"),
                               loess_span = 0.15) {
  smooth_method <- match.arg(smooth_method)

  df_cell <- df_cell %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(time = as.numeric(time))

  y <- df_cell$pct_of_baseline

  y_smooth <- switch(
    smooth_method,
    rollmedian = zoo::rollmedian(y, k = k_smooth, fill = NA, align = "center"),
    rollmean   = zoo::rollmean(y,   k = k_smooth, fill = NA, align = "center"),
    loess = {
      ok <- is.finite(df_cell$time) & is.finite(y)
      if (sum(ok) < 10) {
        rep(NA_real_, length(y))
      } else {
        fit <- stats::loess(y[ok] ~ df_cell$time[ok], span = loess_span, degree = 1)
        out <- rep(NA_real_, length(y))
        out[ok] <- stats::predict(fit, newdata = df_cell$time[ok])
        out
      }
    }
  )

  df_cell %>% dplyr::mutate(pct_smooth = as.numeric(y_smooth))
}

# ---------------------------------------------------------
# Helper: extrema in a window (on pct_smooth by default)
# ---------------------------------------------------------
get_extreme <- function(d, window, which = c("min","max"), col = "pct_smooth") {
  which <- match.arg(which)
  d0 <- d %>% dplyr::filter(time >= window[1], time <= window[2])
  v <- d0[[col]]
  if (nrow(d0) == 0 || all(!is.finite(v))) {
    return(list(time = NA_real_, value = NA_real_))
  }
  idx <- if (which == "min") which.min(v) else which.max(v)
  list(time = d0$time[idx], value = v[idx])
}

# ---------------------------------------------------------
# Focused QC plot: biphasic features (rapid window emphasis)
# - xlim control
# - option to suppress absolute extrema markers
# - option to show ONLY rapid biphasic points (min/max in rapid window)
# ---------------------------------------------------------
plot_biphasic_qc <- function(df_prepped,
                             cell_id,
                             # windows
                             rapid_window = c(0,10),
                             post_window  = c(0,60),
                             # smoothing
                             k_smooth = 11,
                             smooth_method = c("rollmean","rollmedian","loess"),
                             loess_span = 0.15,
                             # calling thresholds
                             min_excursion_pct = 5,
                             # view control
                             x_limits = c(-5, 25),     # <-- adjust this freely
                             show_abs_extrema = FALSE, # <-- turn off abs max/min markers
                             show_rapid_extrema = TRUE,# <-- show rapid min/max markers
                             show_thresholds = TRUE,
                             show_windows = TRUE) {

  smooth_method <- match.arg(smooth_method)

  d <- df_prepped %>%
    dplyr::filter(cell_name == cell_id) %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::arrange(time)

  if (nrow(d) == 0) stop("No rows found for cell_name = ", cell_id)
  if (!"pct_of_baseline" %in% names(d)) stop("df_prepped must include pct_of_baseline")

  d <- add_smoothed_trace(d, k_smooth = k_smooth, smooth_method = smooth_method, loess_span = loess_span)

  # rapid extrema
  rapid_min <- get_extreme(d, rapid_window, "min", col = "pct_smooth")
  rapid_max <- get_extreme(d, rapid_window, "max", col = "pct_smooth")

  # abs extrema (optional)
  abs_min <- get_extreme(d, post_window, "min", col = "pct_smooth")
  abs_max <- get_extreme(d, post_window, "max", col = "pct_smooth")

  # biphasic detection inside rapid window
  rapid_has_inh <- is.finite(rapid_min$value) && rapid_min$value <= (100 - min_excursion_pct)
  rapid_has_exc <- is.finite(rapid_max$value) && rapid_max$value >= (100 + min_excursion_pct)

  biphasic_order <- NA_character_
  if (rapid_has_inh && rapid_has_exc) {
    if (is.finite(rapid_min$time) && is.finite(rapid_max$time)) {
      biphasic_order <- if (rapid_max$time < rapid_min$time) "exc_then_inh" else "inh_then_exc"
    } else {
      biphasic_order <- "ambiguous"
    }
  }

  # colors: RED excitation, BLUE inhibition
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  title_txt <- paste0(
    cell_id,
    " | rapid biphasic: ",
    ifelse(is.na(biphasic_order), "no", biphasic_order)
  )

  subtitle_txt <- paste0(
    "rapid min=", ifelse(is.finite(rapid_min$value), sprintf("%.1f", rapid_min$value), "NA"),
    " @", ifelse(is.finite(rapid_min$time), sprintf("%.1fs", rapid_min$time), "NA"),
    " | rapid max=", ifelse(is.finite(rapid_max$value), sprintf("%.1f", rapid_max$value), "NA"),
    " @", ifelse(is.finite(rapid_max$time), sprintf("%.1fs", rapid_max$time), "NA")
  )

  p <- ggplot(d, aes(x = time, y = pct_of_baseline)) +
    { if (show_windows)
      annotate("rect", xmin = rapid_window[1], xmax = rapid_window[2],
               ymin = -Inf, ymax = Inf, alpha = 0.08)
      else NULL } +
    geom_line(color = "grey60", linewidth = 0.3, alpha = 0.7) +
    geom_line(aes(y = pct_smooth), color = "black", linewidth = 0.8, na.rm = TRUE) +
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5, color = "grey30") +
    { if (show_thresholds) geom_hline(yintercept = 100 + min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_exc) else NULL } +
    { if (show_thresholds) geom_hline(yintercept = 100 - min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_inh) else NULL } +

    # rapid biphasic points (only within rapid window)
    { if (show_rapid_extrema && is.finite(rapid_max$time) && is.finite(rapid_max$value))
      geom_point(data = tibble::tibble(time = rapid_max$time, y = rapid_max$value),
                 aes(x = time, y = y), color = col_exc, size = 2.8)
      else NULL } +
    { if (show_rapid_extrema && is.finite(rapid_min$time) && is.finite(rapid_min$value))
      geom_point(data = tibble::tibble(time = rapid_min$time, y = rapid_min$value),
                 aes(x = time, y = y), color = col_inh, size = 2.8)
      else NULL } +

    # optional absolute extrema markers
    { if (show_abs_extrema && is.finite(abs_max$time) && is.finite(abs_max$value))
      geom_point(data = tibble::tibble(time = abs_max$time, y = abs_max$value),
                 aes(x = time, y = y), shape = 21, fill = col_exc, color = "black", size = 2.8, stroke = 0.4)
      else NULL } +
    { if (show_abs_extrema && is.finite(abs_min$time) && is.finite(abs_min$value))
      geom_point(data = tibble::tibble(time = abs_min$time, y = abs_min$value),
                 aes(x = time, y = y), shape = 21, fill = col_inh, color = "black", size = 2.8, stroke = 0.4)
      else NULL } +

    coord_cartesian(xlim = x_limits) +
    labs(title = title_txt, subtitle = subtitle_txt,
         x = "Time (s)", y = "% of baseline (100 = baseline)") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  p
}

# ---------------------------------------------------------
# Patchwork grid: biphasic exc-first, with x limits + options
# ---------------------------------------------------------
plot_biphasic_excfirst_patchwork <- function(df_prepped,
                                             per_cell2,
                                             subclasses = NULL,
                                             n_max = 24,
                                             ncol = 4,
                                             # pass-through plot args
                                             x_limits = c(-5, 25),
                                             show_abs_extrema = FALSE,
                                             show_rapid_extrema = TRUE,
                                             k_smooth = 11,
                                             smooth_method = "rollmean",
                                             min_excursion_pct = 5,
                                             rapid_window = c(0,10),
                                             post_window = c(0,60)) {

  pick <- per_cell2 %>%
    dplyr::filter(response_type == "biphasic",
                  biphasic_order == "exc_then_inh")

  if (!is.null(subclasses)) {
    pick <- pick %>% dplyr::filter(assigned_subclass %in% subclasses)
  }

  if (nrow(pick) == 0) stop("No biphasic exc_then_inh cells found with current filters.")

  # Optional: strongest-first
  if ("response_peak" %in% names(pick)) {
    pick <- pick %>% dplyr::arrange(dplyr::desc(response_peak))
  }

  pick <- pick %>% dplyr::slice_head(n = n_max)

  plots <- lapply(pick$cell_name, function(cid) {
    plot_biphasic_qc(
      df_prepped = df_prepped,
      cell_id = cid,
      rapid_window = rapid_window,
      post_window = post_window,
      k_smooth = k_smooth,
      smooth_method = smooth_method,
      min_excursion_pct = min_excursion_pct,
      x_limits = x_limits,
      show_abs_extrema = show_abs_extrema,
      show_rapid_extrema = show_rapid_extrema
    ) +
      theme(
        plot.title = element_text(size = 9),
        plot.subtitle = element_text(size = 7),
        axis.title = element_blank()
      )
  })

  wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = paste0("Rapid biphasic (excitation → inhibition) | n = ", nrow(pick)),
      subtitle = paste0("x_limits = [", x_limits[1], ", ", x_limits[2], "]",
                        " | abs extrema shown = ", show_abs_extrema,
                        " | rapid extrema shown = ", show_rapid_extrema),
      theme = theme(plot.title = element_text(face = "bold"))
    )
}

# ---------------------------
# Example usage
# ---------------------------

# A) Single cell, zoomed in, ONLY rapid biphasic points (no abs max/min)
p_one <- plot_biphasic_qc(
  df_prepped, "QF25.26.024.19.07.04",
  x_limits = c(-2, 20),
  show_abs_extrema = FALSE,
  show_rapid_extrema = TRUE,
  k_smooth = 3, smooth_method = "rollmean"
)
print(p_one)

# B) Patchwork: biphasic exc-first, zoom, no abs extrema
p_patch <- plot_biphasic_excfirst_patchwork(
  df_prepped = df_prepped,
  per_cell2  = per_cell2,
  subclasses = c("L23_IT","L5_ET","L5_IT"),
  n_max = 20, ncol = 4,
  x_limits = c(-2, 22),
  show_abs_extrema = FALSE,
  show_rapid_extrema = TRUE,
  k_smooth = 3,
  smooth_method = "rollmean",
  min_excursion_pct = 5
)
print(p_patch)

#ggsave("biphasic_exc_first_patchwork_zoom.png", p_patch, width = 14, height = 10, dpi = 250)
