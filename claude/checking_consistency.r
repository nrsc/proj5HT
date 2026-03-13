## =========================================================
## Puff response QC + bucketed feature extraction + plotting
## Assumes df has columns:
## time, instRate, protocol, percent_change, cell_name, assigned_subclass, ...
##
## Key rules:
## - Filter out protocol == "baseline"
## - Baseline = mean(instRate) in [-5, 0) seconds (per cell)
## - Convert to pct_of_baseline: 100 * instRate / baseline
##   (so 90 = 10% inhibition; 110 = 10% excitation)
## - Buckets:
##   rapid: 0-10 s
##   mid: 10-40 s
##   late_mean: 50-60 s (default)
##   late_slope: 40-60 s (default)
##   abs extrema computed over post_window (default 0-60)
## - Biphasic detection inside rapid window:
##   If both min <= 100-min_excursion_pct and max >= 100+min_excursion_pct occur,
##   label rapid_biphasic and determine order by time of extrema.
##
## Colors (per your latest):
##   RED = excitation
##   BLUE = inhibition
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(zoo)
  library(stringr)
  library(forcats)
  library(scales)
})

## -----------------------------
## 1) Prep raw df
## -----------------------------
# prep_puff_df <- function(df,
#                          baseline_window = c(-5, 0),   # used only if baseline_mode != "auto"
#                          baseline_mode = c("auto", "fixed"),
#                          # windows for auto mode
#                          w_short = c(-5, -1),
#                          w_early = c(-10, -6),
#                          w_long  = c(-10, -1),
#                          # decision thresholds (tune as needed)
#                          z_diff_thr = 1.0,        # |mean_short - mean_early| / pooled_sd
#                          step_hz_thr = 1.0,       # absolute step size (Hz) allowed after slope adjustment
#                          r2_consistent_thr = 0.6, # "consistent slope" indicator
#                          slope_hz_per_s_thr = 0.2,# max allowed |slope| (Hz/s) to still treat as stable baseline
#                          min_pts_per_win = 5) {
#
#   baseline_mode <- match.arg(baseline_mode)
#
#   stopifnot(all(c("cell_name","time","instRate") %in% names(df)))
#
#   df0 <- df %>%
#     dplyr::mutate(
#       time = as.numeric(.data$time),
#       instRate = as.numeric(.data$instRate)
#     )
#
#   # helper: safe stats in a window
#   win_stats <- function(d, w) {
#     d %>%
#       dplyr::filter(.data$time >= w[1], .data$time < w[2]) %>%
#       dplyr::summarise(
#         n = sum(is.finite(.data$instRate)),
#         mean = mean(.data$instRate, na.rm = TRUE),
#         sd = stats::sd(.data$instRate, na.rm = TRUE),
#         cv = dplyr::if_else(is.finite(.data$mean) & .data$mean != 0,
#                             .data$sd / abs(.data$mean),
#                             NA_real_),
#         .groups = "drop"
#       )
#   }
#
#   # helper: linear fit in long window (per cell)
#   long_lm_stats <- function(d, w) {
#     dsub <- d %>% dplyr::filter(.data$time >= w[1], .data$time < w[2]) %>%
#       dplyr::filter(is.finite(.data$time), is.finite(.data$instRate))
#
#     if (nrow(dsub) < 3) {
#       return(dplyr::tibble(slope = NA_real_, r2 = NA_real_))
#     }
#
#     fit <- stats::lm(instRate ~ time, data = dsub)
#     sm  <- summary(fit)
#     dplyr::tibble(
#       slope = unname(stats::coef(fit)[["time"]]),
#       r2 = unname(sm$r.squared)
#     )
#   }
#
#   # --- per-cell diagnostics + baseline choice ---
#   per_cell <- df0 %>%
#     dplyr::group_by(.data$cell_name) %>%
#     dplyr::group_modify(function(d, key) {
#
#       s_short <- win_stats(d, w_short)
#       s_early <- win_stats(d, w_early)
#       s_lm    <- long_lm_stats(d, w_long)
#
#       # pooled SD for a scale-free "separation" indicator
#       pooled_sd <- sqrt((s_short$sd^2 + s_early$sd^2) / 2)
#       z_diff <- dplyr::if_else(is.finite(pooled_sd) & pooled_sd > 0,
#                                abs(s_short$mean - s_early$mean) / pooled_sd,
#                                NA_real_)
#
#       # slope-adjusted "step" between early and short midpoints:
#       # if the difference is just a smooth drift (linear), that's less suspicious.
#       t_mid_short <- mean(w_short)
#       t_mid_early <- mean(w_early)
#       expected_diff_from_slope <- s_lm$slope * (t_mid_short - t_mid_early)
#       step_resid_hz <- (s_short$mean - s_early$mean) - expected_diff_from_slope
#       step_abs_hz <- abs(step_resid_hz)
#
#       # basic validity
#       enough_pts <- (s_short$n >= min_pts_per_win) & (s_early$n >= min_pts_per_win)
#
#       # Decide whether long window is "safe"
#       slope_ok <- is.finite(s_lm$slope) && abs(s_lm$slope) <= slope_hz_per_s_thr
#       r2_ok    <- is.finite(s_lm$r2) && s_lm$r2 >= r2_consistent_thr
#
#       # We FAVOR long window if:
#       # 1) windows cluster (z small) OR step residual small
#       # AND 2) enough points
#       # AND 3) either slope is small OR (slope is consistent and not crazy steep)
#       windows_cluster <- is.finite(z_diff) && (z_diff <= z_diff_thr)
#       step_small      <- is.finite(step_abs_hz) && (step_abs_hz <= step_hz_thr)
#
#       allow_long <- enough_pts &&
#         (windows_cluster || step_small) &&
#         (slope_ok || (r2_ok && abs(s_lm$slope) <= (2 * slope_hz_per_s_thr)))
#
#       # choose baseline window
#       chosen <- if (allow_long) "[-10,-1)" else "[-5,-1)"
#
#       # compute baseline mean from chosen window
#       w_use <- if (allow_long) w_long else w_short
#       s_use <- win_stats(d, w_use)
#
#       reason <- dplyr::case_when(
#         !enough_pts ~ "insufficient_points",
#         allow_long & windows_cluster ~ "favored_long_windows_cluster",
#         allow_long & step_small ~ "favored_long_step_explained_by_slope_or_small",
#         TRUE ~ "fallback_short_possible_rollover_or_separation"
#       )
#
#       dplyr::tibble(
#         baseline_window_used = chosen,
#         baseline_hz = s_use$mean,
#         n_baseline_pts = s_use$n,
#
#         # diagnostics you asked for
#         short_mean = s_short$mean, short_sd = s_short$sd, short_cv = s_short$cv, short_n = s_short$n,
#         early_mean = s_early$mean, early_sd = s_early$sd, early_cv = s_early$cv, early_n = s_early$n,
#
#         mean_diff_short_minus_early = s_short$mean - s_early$mean,
#         pooled_sd = pooled_sd,
#         z_diff = z_diff,
#
#         slope_hz_per_s = s_lm$slope,
#         r2_long = s_lm$r2,
#         step_resid_hz = step_resid_hz,
#         step_abs_hz = step_abs_hz,
#
#         baseline_choice_reason = reason
#       )
#     }) %>%
#     dplyr::ungroup()
#
#   # --- Attach baseline and compute % of baseline ---
#   df1 <- df0 %>%
#     dplyr::left_join(per_cell, by = "cell_name") %>%
#     dplyr::mutate(
#       pct_of_baseline = dplyr::if_else(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
#                                        100 * .data$instRate / .data$baseline_hz,
#                                        NA_real_)
#     )
#
#   if (baseline_mode == "fixed") {
#     # original behavior (overrides auto diagnostics)
#     base_tbl_fixed <- df0 %>%
#       dplyr::filter(.data$time >= baseline_window[1], .data$time < baseline_window[2]) %>%
#       dplyr::group_by(.data$cell_name) %>%
#       dplyr::summarise(
#         baseline_hz = mean(.data$instRate, na.rm = TRUE),
#         n_baseline_pts = sum(is.finite(.data$instRate)),
#         baseline_window_used = paste0("[", baseline_window[1], ",", baseline_window[2], ")"),
#         baseline_choice_reason = "fixed",
#         .groups = "drop"
#       )
#
#     df1 <- df0 %>%
#       dplyr::left_join(base_tbl_fixed, by = "cell_name") %>%
#       dplyr::mutate(
#         pct_of_baseline = dplyr::if_else(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
#                                          100 * .data$instRate / .data$baseline_hz,
#                                          NA_real_)
#       )
#   }
#
#   df1
# }
#
#

prep_puff_df <- function(df,
                         baseline_window = c(-5, 0),   # used only if baseline_mode != "auto"
                         baseline_mode = c("auto", "fixed"),
                         # windows for auto mode
                         w_short = c(-5, -1),
                         w_early = c(-10, -6),
                         w_long  = c(-10, -1),
                         # decision thresholds (tune as needed)
                         z_diff_thr = 1.0,        # |mean_short - mean_early| / pooled_sd
                         step_hz_thr = 1.0,       # absolute step size (Hz) allowed after slope adjustment
                         r2_consistent_thr = 0.6, # "consistent slope" indicator
                         slope_hz_per_s_thr = 0.2,# max allowed |slope| (Hz/s) to still treat as stable baseline
                         min_pts_per_win = 5,

                         # ---- NEW: x_bin anchors ----
                         add_xbin_anchors = TRUE,
                         x_bin_s = 2.5,
                         xlim = c(-10, 50),
                         anchor_fill_instRate = 0,
                         anchor_id_cols = c("cell_name") # or c("cell_name","protocol") if you want
) {

  baseline_mode <- match.arg(baseline_mode)

  stopifnot(all(c("cell_name","time","instRate") %in% names(df)))

  df0 <- df %>%
    dplyr::mutate(
      time = as.numeric(.data$time),
      instRate = as.numeric(.data$instRate)
    )

  # helper: safe stats in a window
  win_stats <- function(d, w) {
    d %>%
      dplyr::filter(.data$time >= w[1], .data$time < w[2]) %>%
      dplyr::summarise(
        n = sum(is.finite(.data$instRate)),
        mean = mean(.data$instRate, na.rm = TRUE),
        sd = stats::sd(.data$instRate, na.rm = TRUE),
        cv = dplyr::if_else(is.finite(.data$mean) & .data$mean != 0,
                            .data$sd / abs(.data$mean),
                            NA_real_),
        .groups = "drop"
      )
  }

  # helper: linear fit in long window (per cell)
  long_lm_stats <- function(d, w) {
    dsub <- d %>%
      dplyr::filter(.data$time >= w[1], .data$time < w[2]) %>%
      dplyr::filter(is.finite(.data$time), is.finite(.data$instRate))

    if (nrow(dsub) < 3) {
      return(dplyr::tibble(slope = NA_real_, r2 = NA_real_))
    }

    fit <- stats::lm(instRate ~ time, data = dsub)
    sm  <- summary(fit)
    dplyr::tibble(
      slope = unname(stats::coef(fit)[["time"]]),
      r2 = unname(sm$r.squared)
    )
  }

  # --- per-cell diagnostics + baseline choice ---
  per_cell <- df0 %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(function(d, key) {

      s_short <- win_stats(d, w_short)
      s_early <- win_stats(d, w_early)
      s_lm    <- long_lm_stats(d, w_long)

      pooled_sd <- sqrt((s_short$sd^2 + s_early$sd^2) / 2)
      z_diff <- dplyr::if_else(is.finite(pooled_sd) & pooled_sd > 0,
                               abs(s_short$mean - s_early$mean) / pooled_sd,
                               NA_real_)

      t_mid_short <- mean(w_short)
      t_mid_early <- mean(w_early)
      expected_diff_from_slope <- s_lm$slope * (t_mid_short - t_mid_early)
      step_resid_hz <- (s_short$mean - s_early$mean) - expected_diff_from_slope
      step_abs_hz <- abs(step_resid_hz)

      enough_pts <- (s_short$n >= min_pts_per_win) & (s_early$n >= min_pts_per_win)

      slope_ok <- is.finite(s_lm$slope) && abs(s_lm$slope) <= slope_hz_per_s_thr
      r2_ok    <- is.finite(s_lm$r2) && s_lm$r2 >= r2_consistent_thr

      windows_cluster <- is.finite(z_diff) && (z_diff <= z_diff_thr)
      step_small      <- is.finite(step_abs_hz) && (step_abs_hz <= step_hz_thr)

      allow_long <- enough_pts &&
        (windows_cluster || step_small) &&
        (slope_ok || (r2_ok && abs(s_lm$slope) <= (2 * slope_hz_per_s_thr)))

      chosen <- if (allow_long) "[-10,-1)" else "[-5,-1)"
      w_use  <- if (allow_long) w_long else w_short
      s_use  <- win_stats(d, w_use)

      reason <- dplyr::case_when(
        !enough_pts ~ "insufficient_points",
        allow_long & windows_cluster ~ "favored_long_windows_cluster",
        allow_long & step_small ~ "favored_long_step_explained_by_slope_or_small",
        TRUE ~ "fallback_short_possible_rollover_or_separation"
      )

      dplyr::tibble(
        baseline_window_used = chosen,
        baseline_hz = s_use$mean,
        n_baseline_pts = s_use$n,

        short_mean = s_short$mean, short_sd = s_short$sd, short_cv = s_short$cv, short_n = s_short$n,
        early_mean = s_early$mean, early_sd = s_early$sd, early_cv = s_early$cv, early_n = s_early$n,

        mean_diff_short_minus_early = s_short$mean - s_early$mean,
        pooled_sd = pooled_sd,
        z_diff = z_diff,

        slope_hz_per_s = s_lm$slope,
        r2_long = s_lm$r2,
        step_resid_hz = step_resid_hz,
        step_abs_hz = step_abs_hz,

        baseline_choice_reason = reason
      )
    }) %>%
    dplyr::ungroup()

  # --- Attach baseline and compute % of baseline ---
  df1 <- df0 %>%
    dplyr::left_join(per_cell, by = "cell_name") %>%
    dplyr::mutate(
      pct_of_baseline = dplyr::if_else(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
                                       100 * .data$instRate / .data$baseline_hz,
                                       NA_real_)
    )

  if (baseline_mode == "fixed") {
    base_tbl_fixed <- df0 %>%
      dplyr::filter(.data$time >= baseline_window[1], .data$time < baseline_window[2]) %>%
      dplyr::group_by(.data$cell_name) %>%
      dplyr::summarise(
        baseline_hz = mean(.data$instRate, na.rm = TRUE),
        n_baseline_pts = sum(is.finite(.data$instRate)),
        baseline_window_used = paste0("[", baseline_window[1], ",", baseline_window[2], ")"),
        baseline_choice_reason = "fixed",
        .groups = "drop"
      )

    df1 <- df0 %>%
      dplyr::left_join(base_tbl_fixed, by = "cell_name") %>%
      dplyr::mutate(
        pct_of_baseline = dplyr::if_else(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
                                         100 * .data$instRate / .data$baseline_hz,
                                         NA_real_)
      )
  }

  # ---- NEW: add x_bin + anchor rows here (ideal place) ----
  if (isTRUE(add_xbin_anchors)) {
    stopifnot(all(anchor_id_cols %in% names(df1)))

    x_grid <- seq(from = xlim[1], to = xlim[2], by = x_bin_s)

    df1 <- df1 %>%
      dplyr::mutate(
        # choose your binning rule; this matches "bin start" convention
        x_bin = floor(.data$time / x_bin_s) * x_bin_s
        # if you really want your +5 shift, use:
        # x_bin = floor(.data$time / x_bin_s) * x_bin_s + 5
      ) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(anchor_id_cols))) %>%
      tidyr::complete(x_bin = x_grid) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        is_anchor = is.na(.data$time) & is.na(.data$instRate),
        instRate = dplyr::coalesce(.data$instRate, anchor_fill_instRate),

        # keep baseline_hz carried through; recompute pct_of_baseline after fill
        pct_of_baseline = dplyr::if_else(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
                                         100 * .data$instRate / .data$baseline_hz,
                                         NA_real_)
      ) %>%
      # fill metadata columns for the new rows so facets/grouping keep working
      dplyr::group_by(dplyr::across(dplyr::all_of(anchor_id_cols))) %>%
      tidyr::fill(dplyr::everything(), .direction = "downup") %>%
      dplyr::ungroup()
  }

  df1
}




## -----------------------------
## 2) Smoothing helper
## -----------------------------
add_smoothed_trace <- function(df_cell,
                               k_smooth = 11,
                               smooth_method = c("rollmedian","rollmean","loess"),
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

## -----------------------------
## Cell Binning
# --- Adaptive bin widths tuned to your baseline_hz range (~2–13 Hz) ---
choose_bins_for_cell <- function(dt, baseline_hz) {

  # RAPID: preserve timing for biphasic order
  bin_rapid <- 0.5

  # POST/MID/LATE: rate-aware, stable but not too chunky
  target_spikes <- 4
  bin_post <- if (is.finite(baseline_hz) && baseline_hz > 0) target_spikes / baseline_hz else 1.0
  bin_post <- max(0.5, min(1.5, bin_post))  # clamp

  list(
    rapid = max(bin_rapid, dt),
    post  = max(bin_post,  dt),
    mid   = max(bin_post,  dt),
    late  = max(bin_post,  dt)
  )
}


## -----------------------------
## inhibition rescue
inh_area_binned <- function(dxy, thr, window = c(0,10), use_col = c("y_s","y")) {
  use_col <- match.arg(use_col)
  if (nrow(dxy) < 2) return(0)

  d <- dxy %>% dplyr::filter(time >= window[1], time <= window[2])
  y <- d[[use_col]]
  ok <- is.finite(d$time) & is.finite(y)
  d <- d[ok, , drop = FALSE]
  if (nrow(d) < 2) return(0)

  # area of excursion below (100-thr), in "%*s" units
  deficit <- pmax(0, (100 - thr) - y)
  dt <- c(diff(d$time), stats::median(diff(d$time), na.rm=TRUE))
  dt[!is.finite(dt)] <- 0
  sum(deficit * dt, na.rm = TRUE)
}

# --- Bin a trace in time and optionally smooth within the binned series ---
# --- Bin a trace in time and optionally smooth within the binned series ---
# Key fixes:
# 1) Use floor()-based binning anchored at window[1] (prevents duplicate bins)
# 2) Use dplyr::reframe() (silences dplyr 1.1+ summarise size warnings)
# 3) Robust dt/k computation + edge fill so early/late bins remain usable
bin_and_smooth <- function(df_cell, window, time_bin_s,
                           y_col = "pct_of_baseline",
                           bin_fun = c("median","mean"),
                           smooth_s = 0.5,              # smoothing in SECONDS on the binned series
                           smooth_fun = c("rollmean","rollmedian"),
                           digits_time = 6) {           # optional: rounding for stable printing/joins

  bin_fun <- match.arg(bin_fun)
  smooth_fun <- match.arg(smooth_fun)

  stopifnot(length(window) == 2, is.numeric(window), is.finite(window[1]), is.finite(window[2]))
  stopifnot(is.numeric(time_bin_s), length(time_bin_s) == 1, is.finite(time_bin_s), time_bin_s > 0)

  if (!all(c("time", y_col) %in% names(df_cell))) {
    stop("bin_and_smooth: df_cell must contain columns: time and ", y_col)
  }

  d0 <- df_cell %>%
    dplyr::mutate(time = as.numeric(.data$time)) %>%
    dplyr::filter(.data$time >= window[1], .data$time <= window[2]) %>%
    dplyr::filter(is.finite(.data$time), is.finite(.data[[y_col]])) %>%
    dplyr::arrange(.data$time)

  if (nrow(d0) == 0) {
    return(tibble::tibble(
      time = numeric(0), y = numeric(0), y_s = numeric(0), n = integer(0),
      k_smooth_bins = integer(0), dt_bin = numeric(0)
    ))
  }

  # ---- robust, duplicate-free binning ----
  # anchor bins at window start (t0) and use floor -> deterministic
  t0 <- window[1]
  d_bin <- d0 %>%
    dplyr::mutate(
      time_bin = t0 + floor((.data$time - t0) / time_bin_s) * time_bin_s,
      time_bin = round(.data$time_bin, digits_time)  # optional rounding for stability
    ) %>%
    dplyr::group_by(.data$time_bin) %>%
    dplyr::reframe(
      time = .data$time_bin[1],
      y = if (bin_fun == "median") stats::median(.data[[y_col]], na.rm = TRUE)
      else mean(.data[[y_col]], na.rm = TRUE),
      n = dplyr::n()
    ) %>%
    dplyr::arrange(.data$time)

  # ---- smoothing on the BINNED series in seconds (not points) ----
  dtb <- stats::median(diff(d_bin$time), na.rm = TRUE)
  if (!is.finite(dtb) || dtb <= 0) dtb <- time_bin_s

  k <- max(1L, as.integer(round(smooth_s / dtb)))
  if (k %% 2L == 0L) k <- k + 1L  # prefer odd window for symmetry

  y_s <- if (k <= 1L) {
    d_bin$y
  } else if (smooth_fun == "rollmean") {
    as.numeric(zoo::rollmean(d_bin$y, k = k, fill = NA, align = "center"))
  } else {
    as.numeric(zoo::rollmedian(d_bin$y, k = k, fill = NA, align = "center"))
  }

  # edge fill: keep endpoints usable
  y_s[!is.finite(y_s)] <- d_bin$y[!is.finite(y_s)]

  d_bin %>%
    dplyr::mutate(
      y_s = y_s,
      k_smooth_bins = k,
      dt_bin = dtb
    )
}


# --- Persistence on a (time, y_s) binned+smoothed series ---
has_persist_binned <- function(dxy, direction = c("exc","inh"),
                               thr = 10, min_persist_s = 0.5,
                               use_col = c("y_s","y")) {
  direction <- match.arg(direction)
  use_col <- match.arg(use_col)

  if (nrow(dxy) < 2) return(FALSE)

  # choose series; if y_s is too sparse, fall back to y
  y <- dxy[[use_col]]
  ok <- is.finite(dxy$time) & is.finite(y)

  if (sum(ok) < 2 && use_col == "y_s" && "y" %in% names(dxy)) {
    y <- dxy$y
    ok <- is.finite(dxy$time) & is.finite(y)
  }

  d2 <- dxy[ok, , drop = FALSE]
  if (nrow(d2) < 2) return(FALSE)

  dtb <- stats::median(diff(d2$time), na.rm = TRUE)
  if (!is.finite(dtb) || dtb <= 0) dtb <- 0

  flag <- if (direction == "exc") y[ok] >= (100 + thr) else y[ok] <= (100 - thr)
  if (!any(flag)) return(FALSE)

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  true_runs <- which(r$values)

  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]
    dur <- (d2$time[i1] - d2$time[i0]) + dtb
    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }

  FALSE
}

# --- Extremum on binned+smoothed series ---
get_extreme_binned <- function(dxy, which = c("min","max"), col = c("y_s","y")) {
  which <- match.arg(which)
  col <- match.arg(col)

  if (nrow(dxy) == 0 || all(!is.finite(dxy[[col]]))) return(list(time = NA_real_, value = NA_real_))
  v <- dxy[[col]]
  idx <- if (which == "min") which.min(v) else which.max(v)
  list(time = dxy$time[idx], value = v[idx])
}



## -----------------------------
## 2b) Persistence helper
## -----------------------------
# returns TRUE if there exists a continuous run (in time) beyond threshold
# of at least min_persist_s seconds inside the window.
has_persistent_excursion <- function(d, window, direction = c("exc","inh"),
                                     thr = 10,
                                     col = "pct_smooth",
                                     min_persist_s = 0.5,
                                     dt_fallback = NULL) {
  direction <- match.arg(direction)

  d0 <- d %>%
    dplyr::filter(time >= window[1], time <= window[2]) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(time))

  if (nrow(d0) < 2) return(FALSE)

  # estimate sampling interval (robust)
  dt <- stats::median(diff(d0$time), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) dt <- dt_fallback
  if (!is.finite(dt) || dt <= 0) dt <- 0  # last resort

  flag <- if (direction == "exc") d0[[col]] >= (100 + thr) else d0[[col]] <= (100 - thr)
  if (!any(flag)) return(FALSE)

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  true_runs <- which(r$values)
  if (length(true_runs) == 0) return(FALSE)

  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]

    # IMPORTANT: add ~1 sample to duration so single-point events aren’t dur=0
    dur <- (d0$time[i1] - d0$time[i0]) + dt

    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }

  FALSE
}

# longest continuous run duration (seconds) where instRate <= zero_thr
max_zero_run_s <- function(time, instRate, window = c(0, 10), zero_thr = 0, dt_fallback = 0.1) {
  ok <- is.finite(time) & is.finite(instRate)
  time <- as.numeric(time[ok])
  instRate <- as.numeric(instRate[ok])

  keep <- time >= window[1] & time <= window[2]
  time <- time[keep]; instRate <- instRate[keep]
  if (length(time) < 2) return(0)

  o <- order(time)
  time <- time[o]; instRate <- instRate[o]

  dt <- stats::median(diff(time), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) dt <- dt_fallback

  flag <- instRate <= zero_thr
  if (!any(flag)) return(0)

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  true_runs <- which(r$values)

  max_dur <- 0
  for (k in true_runs) {
    i0 <- starts[k]; i1 <- ends[k]
    dur <- (time[i1] - time[i0]) + dt
    if (is.finite(dur)) max_dur <- max(max_dur, dur)
  }

  max_dur
}


## -----------------------------
## 3) Feature extraction per cell
## -----------------------------
extract_bucket_features_one <- function(df_cell,
                                        rapid_window = c(0, 10),
                                        mid_window   = c(10, 40),
                                        late_slope_window = c(40, 60),
                                        late_mean_window  = c(50, 60),
                                        post_window = c(0, 60),
                                        min_excursion_pct = 5,
                                        k_smooth = 11,
                                        smooth_method = c("rollmedian","rollmean","loess"),
                                        loess_span = 0.15,
                                        min_persist_s = 0.5,
                                        zero_inh_s = 3,
                                        zero_hz_thr = 0
) {

  smooth_method <- match.arg(smooth_method)

  df_cell <- df_cell %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::arrange(time)

  if (!"pct_of_baseline" %in% names(df_cell)) {
    stop("extract_bucket_features_one: missing pct_of_baseline. Run prep_puff_df() first.")
  }

  dt <- stats::median(diff(sort(unique(df_cell$time))), na.rm = TRUE)
  if (!is.finite(dt) || dt <= 0) dt <- 0.25

  baseline_hz_val <- if ("baseline_hz" %in% names(df_cell)) df_cell$baseline_hz[1] else NA_real_
  rapid_zero_run_s <- max_zero_run_s(
    time = df_cell$time,
    instRate = df_cell$instRate,      # uses RAW spiking
    window = rapid_window,
    zero_thr = zero_hz_thr,
    dt_fallback = max(0.05, dt)
  )

  rapid_zero_inh <- is.finite(rapid_zero_run_s) && (rapid_zero_run_s >= zero_inh_s)


  bins <- choose_bins_for_cell(dt, baseline_hz_val)

  d_rapid <- bin_and_smooth(df_cell, rapid_window, time_bin_s = bins$rapid,
                            y_col = "pct_of_baseline", bin_fun = "median",
                            smooth_s = 0.5, smooth_fun = "rollmean")

  d_post  <- bin_and_smooth(df_cell, post_window,  time_bin_s = bins$post,
                            y_col = "pct_of_baseline", bin_fun = "median",
                            smooth_s = 0.75, smooth_fun = "rollmean")

  d_mid   <- bin_and_smooth(df_cell, mid_window,   time_bin_s = bins$mid,
                            y_col = "pct_of_baseline", bin_fun = "median",
                            smooth_s = 1.0, smooth_fun = "rollmean")

  d_lateS <- bin_and_smooth(df_cell, late_slope_window, time_bin_s = bins$late,
                            y_col = "pct_of_baseline", bin_fun = "median",
                            smooth_s = 1.0, smooth_fun = "rollmean")

  d_lateM <- bin_and_smooth(df_cell, late_mean_window,  time_bin_s = bins$late,
                            y_col = "pct_of_baseline", bin_fun = "median",
                            smooth_s = 1.0, smooth_fun = "rollmean")

  rapid_min <- get_extreme_binned(d_rapid, "min", col = "y_s")
  rapid_max <- get_extreme_binned(d_rapid, "max", col = "y_s")
  abs_min   <- get_extreme_binned(d_post,  "min", col = "y_s")
  abs_max   <- get_extreme_binned(d_post,  "max", col = "y_s")

  rapid_has_exc_persist <- has_persist_binned(d_rapid, "exc",
                                              thr = min_excursion_pct,
                                              min_persist_s = min_persist_s,
                                              use_col = "y_s")

  rapid_has_inh_persist <- has_persist_binned(d_rapid, "inh",
                                              thr = min_excursion_pct,
                                              min_persist_s = min_persist_s,
                                              use_col = "y_s")

  abs_has_exc_persist <- has_persist_binned(d_post, "exc",
                                            thr = min_excursion_pct,
                                            min_persist_s = min_persist_s,
                                            use_col = "y_s")

  abs_has_inh_persist <- has_persist_binned(d_post, "inh",
                                            thr = min_excursion_pct,
                                            min_persist_s = min_persist_s,
                                            use_col = "y_s")

  max_is_real <- is.finite(abs_max$value) &&
    abs_max$value >= (100 + min_excursion_pct) &&
    abs_has_exc_persist

  abs_min_ok_time <- is.finite(abs_min$time) && abs_min$time > max(0.1, bins$post)

  min_is_real <- abs_min_ok_time &&
    is.finite(abs_min$value) &&
    abs_min$value <= (100 - min_excursion_pct) &&
    abs_has_inh_persist

  abs_max_use <- if (max_is_real) abs_max$value else NA_real_
  abs_min_use <- if (min_is_real) abs_min$value else NA_real_

  exc_mag <- if (is.finite(abs_max_use)) abs_max_use - 100 else NA_real_
  inh_mag <- if (is.finite(abs_min_use)) 100 - abs_min_use else NA_real_

  direction <- dplyr::case_when(
    is.finite(exc_mag) & is.finite(inh_mag) ~ ifelse(exc_mag >= inh_mag, "excitation", "inhibition"),
    is.finite(exc_mag) & !is.finite(inh_mag) ~ "excitation",
    !is.finite(exc_mag) & is.finite(inh_mag) ~ "inhibition",
    TRUE ~ NA_character_
  )

  response_peak <- dplyr::case_when(
    direction == "excitation" ~ exc_mag,
    direction == "inhibition" ~ inh_mag,
    TRUE ~ NA_real_
  )

  rapid_pattern <- dplyr::case_when(
    rapid_has_exc_persist & rapid_has_inh_persist ~
      ifelse(rapid_max$time < rapid_min$time, "rapid_biphasic_exc_to_inh", "rapid_biphasic_inh_to_exc"),
    rapid_has_exc_persist ~ "rapid_excitation_only",
    rapid_has_inh_persist ~ "rapid_inhibition_only",
    TRUE ~ "rapid_no_clear_change"
  )

  # ---- MID / LATE summary from SMOOTHED (existing) ----
  mid_mean  <- if (nrow(d_mid)   > 0) mean(d_mid$y_s,   na.rm = TRUE) else NA_real_
  late_mean <- if (nrow(d_lateM) > 0) mean(d_lateM$y_s, na.rm = TRUE) else NA_real_

  # ---- NEW: MID / LATE summary from RAW-BINNED (semi-smooth support) ----
  mid_mean_raw  <- if (nrow(d_mid)   > 0) mean(d_mid$y,   na.rm = TRUE) else NA_real_
  late_mean_raw <- if (nrow(d_lateM) > 0) mean(d_lateM$y, na.rm = TRUE) else NA_real_
  mid_max_raw   <- if (nrow(d_mid)   > 0) max(d_mid$y,    na.rm = TRUE) else NA_real_

  late_slope <- {
    dlate <- d_lateS %>% dplyr::filter(is.finite(time), is.finite(y_s))
    if (nrow(dlate) >= 5) coef(stats::lm(y_s ~ time, data = dlate))[["time"]] else NA_real_
  }

  rapid_inh_area <- inh_area_binned(d_rapid, thr = min_excursion_pct, window = rapid_window, use_col = "y_s")

  tibble::tibble(
    cell_name = df_cell$cell_name[1],
    assigned_subclass = if ("assigned_subclass" %in% names(df_cell)) df_cell$assigned_subclass[1] else NA_character_,
    baseline_hz = baseline_hz_val,

    rapid_min_time = rapid_min$time,
    rapid_min_pct  = rapid_min$value,
    rapid_max_time = rapid_max$time,
    rapid_max_pct  = rapid_max$value,

    abs_min_time = abs_min$time,
    abs_min_pct  = abs_min$value,
    abs_max_time = abs_max$time,
    abs_max_pct  = abs_max$value,

    min_is_real = min_is_real,
    max_is_real = max_is_real,
    exc_mag = exc_mag,
    inh_mag = inh_mag,
    direction = direction,
    response_peak = response_peak,

    rapid_pattern = rapid_pattern,
    rapid_has_exc_persist = rapid_has_exc_persist,
    rapid_has_inh_persist = rapid_has_inh_persist,

    mid_mean_pct = mid_mean,
    late_mean_pct = late_mean,
    late_slope_pct_per_s = late_slope,

    # NEW outputs
    mid_mean_pct_raw  = mid_mean_raw,
    late_mean_pct_raw = late_mean_raw,
    mid_max_pct_raw   = mid_max_raw,

    k_smooth = k_smooth,
    smooth_method = smooth_method,
    min_persist_s = min_persist_s,

    bin_rapid_s = bins$rapid,
    bin_post_s  = bins$post,
    bin_mid_s   = bins$mid,
    bin_late_s  = bins$late,
    rapid_bin_k = if (nrow(d_rapid)>0) d_rapid$k_smooth_bins[1] else NA_integer_,
    post_bin_k  = if (nrow(d_post)>0)  d_post$k_smooth_bins[1]  else NA_integer_,

    rapid_min_pct_bin  = rapid_min$value,
    rapid_max_pct_bin  = rapid_max$value,
    abs_min_pct_bin    = abs_min$value,
    abs_max_pct_bin    = abs_max$value,

    abs_has_inh_persist = abs_has_inh_persist,
    abs_has_exc_persist = abs_has_exc_persist,
    rapid_inh_area = rapid_inh_area,
    rapid_zero_run_s = rapid_zero_run_s,
    rapid_zero_inh   = rapid_zero_inh

  )
}



extract_bucket_features <- function(df_prepped,
                                    rapid_window = c(0, 10),
                                    mid_window   = c(10, 40),
                                    late_slope_window = c(40, 60),
                                    late_mean_window  = c(50, 60),
                                    post_window = c(0, 60),
                                    min_excursion_pct = 5,
                                    k_smooth = 11,
                                    smooth_method = c("rollmedian","rollmean","loess"),
                                    loess_span = 0.15,
                                    min_persist_s = 0.5) {

  smooth_method <- match.arg(smooth_method)

  df_prepped %>%
    dplyr::group_by(cell_name) %>%
    dplyr::group_modify(~extract_bucket_features_one(
      .x,
      rapid_window = rapid_window,
      mid_window = mid_window,
      late_slope_window = late_slope_window,
      late_mean_window = late_mean_window,
      post_window = post_window,
      min_excursion_pct = min_excursion_pct,
      k_smooth = k_smooth,
      smooth_method = smooth_method,
      loess_span = loess_span,
      min_persist_s = min_persist_s
    )) %>%
    dplyr::ungroup()
}

#extract_bucket_features <- function(df_prepped, ...) {
  df_prepped %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(~extract_bucket_features_one_simple(.x, ...), .keep = TRUE) %>%
    dplyr::ungroup()
}



## -----------------------------
## plot dataframes for modalitites
## -----------------------------

make_plot_df <- function(df_prepped,
                         per_cell,
                         xlim = c(-5, 56),
                         k_smooth = 3,
                         drop_no_change = TRUE) {

  per_cell0 <- per_cell %>%
    dplyr::mutate(response_type = as.character(response_type))

  # ensure biphasic_order exists (avoid if_else on scalar condition)
  if (!"biphasic_order" %in% names(per_cell0)) {
    per_cell0 <- per_cell0 %>% dplyr::mutate(biphasic_order = NA_character_)
  } else {
    per_cell0 <- per_cell0 %>% dplyr::mutate(biphasic_order = as.character(biphasic_order))
  }

  keep_tbl <- per_cell0 %>%
    dplyr::filter(!is.na(response_type)) %>%
    { if (drop_no_change) dplyr::filter(., response_type != "no_change") else . } %>%
    dplyr::mutate(
      biphasic_order = as.character(biphasic_order)
    ) %>%
    # drop unknown biphasic order ONLY for biphasic rows
    dplyr::filter(
      !(response_type == "biphasic" & is.na(biphasic_order))
    ) %>%
    dplyr::select(cell_name, response_type, assigned_subclass, biphasic_order) %>%
    dplyr::distinct()


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

  if (nrow(d) == 0) stop("make_plot_df: no rows after join/xlim filters.")

  d
}

make_plot_df_ei <- function(df_prepped, per_cell2,
                            xlim = c(-5, 56),
                            k_smooth = 3,
                            time_bin_s = NULL,          # if NULL, auto from dt_med_overall
                            min_cells_per_time = 3) {

  # pick a reasonable bin automatically from your data resolution
  if (is.null(time_bin_s)) {
    dt_est <- df_prepped %>%
      dplyr::group_by(cell_name) %>%
      dplyr::summarise(dt = stats::median(diff(sort(unique(time)))), .groups="drop") %>%
      dplyr::summarise(dt = stats::median(dt, na.rm=TRUE)) %>%
      dplyr::pull(dt)
    time_bin_s <- max(0.1, dt_est)  # floor at 100 ms
  }

  # join labels onto the timepoint df
  d0 <- df_prepped %>%
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%  # avoid duplicate column from df_prepped
    dplyr::inner_join(
      per_cell2 %>% dplyr::select(cell_name, response_type, assigned_subclass),
      by = "cell_name"
    ) %>%
    dplyr::filter(
      response_type %in% c("excitation", "inhibition"),
      time >= xlim[1], time <= xlim[2],
      is.finite(pct_of_baseline)
    ) %>%
    dplyr::mutate(time = as.numeric(time))

  # per-cell smoothing (grey traces)
  d_traces <- d0 %>%
    dplyr::group_by(response_type, assigned_subclass, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  # mean per group per time (BLACK mean line) — IMPORTANT: not duplicated per cell
  d_mean <- d_traces %>%
    dplyr::mutate(time_bin = round(time / time_bin_s) * time_bin_s) %>%
    dplyr::group_by(response_type, assigned_subclass, time_bin) %>%
    dplyr::summarise(
      n_cells = dplyr::n_distinct(cell_name[is.finite(pct_smooth)]),
      mean_pct = mean(pct_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(mean_pct), n_cells >= min_cells_per_time) %>%
    dplyr::rename(time = time_bin)

  list(traces = d_traces, mean = d_mean, time_bin_s = time_bin_s)
}



make_mean_df <- function(d,
                         mean_smooth = c("spline","rollmean","none"),
                         spline_df = 20,
                         mean_roll_k = 11) {

  mean_smooth <- match.arg(mean_smooth)

  mean_df <- d %>%
    dplyr::filter(is.finite(pct_smooth), is.finite(time)) %>%
    dplyr::group_by(panel, assigned_subclass, time) %>%
    dplyr::summarise(
      n_cells = dplyr::n_distinct(cell_name),
      mean_y = mean(pct_smooth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(panel, assigned_subclass) %>%
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

  mean_df
}

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

## -----------------------------
## 4) QC plot for a single cell
## -----------------------------
plot_bucket_qc <- function(df_prepped,
                           cell_id,
                           rapid_window = c(0,10),
                           mid_window = c(10,40),
                           late_slope_window = c(40,60),
                           late_mean_window = c(50,60),
                           post_window = c(0,60),
                           k_smooth = 11,
                           smooth_method = c("rollmedian","rollmean","loess"),
                           loess_span = 0.15,
                           min_excursion_pct = 5,
                           show_instRate = FALSE) {

  smooth_method <- match.arg(smooth_method)

  d <- df_prepped %>%
    dplyr::filter(cell_name == cell_id) %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::arrange(time)

  if (nrow(d) == 0) stop("No rows found for cell_name = ", cell_id)
  if (!"pct_of_baseline" %in% names(d)) stop("df_prepped must include pct_of_baseline (run prep_puff_df).")

  d <- add_smoothed_trace(d, k_smooth = k_smooth, smooth_method = smooth_method, loess_span = loess_span)

  feats <- extract_bucket_features_one(
    d,
    rapid_window = rapid_window,
    mid_window = mid_window,
    late_slope_window = late_slope_window,
    late_mean_window = late_mean_window,
    post_window = post_window,
    min_excursion_pct = min_excursion_pct,
    k_smooth = k_smooth,
    smooth_method = smooth_method,
    loess_span = loess_span
  )

  # Late fit segment for plotting
  dlate <- d %>%
    dplyr::filter(time >= late_slope_window[1], time <= late_slope_window[2]) %>%
    dplyr::filter(is.finite(time), is.finite(pct_smooth))

  fit_line <- NULL
  if (nrow(dlate) >= 5) {
    fit <- stats::lm(pct_smooth ~ time, data = dlate)
    fit_line <- tibble::tibble(
      time = seq(min(dlate$time), max(dlate$time), length.out = 50)
    ) %>%
      dplyr::mutate(pct_fit = stats::predict(fit, newdata = .))
  }

  # Marker colors: RED excitation, BLUE inhibition
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  # Choose point colors for abs extrema
  abs_max_col <- col_exc
  abs_min_col <- col_inh

  # Build title
  subclass <- feats$assigned_subclass
  basehz <- feats$baseline_hz
  rp <- feats$rapid_pattern
  dir <- feats$direction
  peak <- feats$response_peak
  title_txt <- paste0(cell_id,
                      if (!is.na(subclass)) paste0(" | ", subclass) else "",
                      if (!is.na(dir)) paste0(" | ", dir) else "",
                      if (!is.na(peak)) paste0(" | peak = ", sprintf("%.1f", peak), "%") else "")

  subtitle_txt <- paste0(
    "Baseline = ", ifelse(is.finite(basehz), sprintf("%.2f Hz", basehz), "NA"),
    " | rapid = ", rp,
    " | abs_min=", ifelse(is.finite(feats$abs_min_pct), sprintf("%.1f%%", feats$abs_min_pct), "NA/flagged"),
    " | abs_max=", ifelse(is.finite(feats$abs_max_pct), sprintf("%.1f%%", feats$abs_max_pct), "NA/flagged")
  )

  # Plot (pct baseline)
  p <- ggplot(d, aes(x = time, y = pct_of_baseline)) +
    # shaded windows
    annotate("rect", xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf, alpha = 0.08) +
    annotate("rect", xmin = mid_window[1], xmax = mid_window[2],
             ymin = -Inf, ymax = Inf, alpha = 0.05) +
    annotate("rect", xmin = late_slope_window[1], xmax = late_slope_window[2],
             ymin = -Inf, ymax = Inf, alpha = 0.06) +
    # raw trace
    geom_line(color = "grey60", linewidth = 0.3, alpha = 0.7) +
    # smoothed trace
    geom_line(aes(y = pct_smooth), color = "black", linewidth = 0.7, na.rm = TRUE) +
    # baseline line at 100
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5, color = "grey30") +
    # thresholds for "real" excursion
    geom_hline(yintercept = 100 + min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_exc) +
    geom_hline(yintercept = 100 - min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_inh) +
    # absolute extrema markers (only if "real")
    { if (isTRUE(feats$max_is_real) && is.finite(feats$abs_max_time) && is.finite(feats$abs_max_pct))
      geom_point(data = tibble::tibble(time = feats$abs_max_time, y = feats$abs_max_pct),
                 aes(x = time, y = y), color = abs_max_col, size = 2.6)
      else NULL } +
    { if (isTRUE(feats$min_is_real) && is.finite(feats$abs_min_time) && is.finite(feats$abs_min_pct))
      geom_point(data = tibble::tibble(time = feats$abs_min_time, y = feats$abs_min_pct),
                 aes(x = time, y = y), color = abs_min_col, size = 2.6)
      else NULL } +
    # late fit line
    { if (!is.null(fit_line))
      geom_line(data = fit_line, aes(x = time, y = pct_fit), linewidth = 1.0)
      else NULL } +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "Time (s)",
      y = "% of baseline (100 = baseline)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  p
}

## -----------------------------
## 5) Loop helper: save many QC plots
## -----------------------------
plot_bucket_qc_many <- function(df_prepped,
                                cell_ids,
                                out_dir = "qc_bucket_plots",
                                width = 7, height = 4.5, dpi = 200,
                                ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  for (cid in cell_ids) {
    p <- plot_bucket_qc(df_prepped, cid, ...)
    fn <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", cid), ".png"))
    ggsave(fn, p, width = width, height = height, dpi = dpi)
  }

  message("Saved ", length(cell_ids), " QC plots to: ", out_dir)
}





# =========================================================
# Integrate biphasic feature into per-cell table
# Uses rapid_* and abs_* columns (NOT early_*).
# =========================================================
add_biphasic_features_consistent <- function(
    per_cell,
    min_excursion_pct = 15,
    rescue_biphasic   = TRUE,
    rescue_inh_window = c(-1, 10),
    rescue_inh_extra  = 0,
    rescue_area_thr   = 1,
    min_sep_s         = 0.25,
    rel_inh_to_exc    = 0.35,
    deep_dip_pct      = 70,
    rebound_mid_thr   = 5,
    raw_rebound_mid_thr = 2,
    raw_use = c("mean","max"),
    still_inh_mid_thr = 10
) {
  raw_use <- match.arg(raw_use)

  has_cols <- function(...) all(c(...) %in% names(per_cell))
  need <- c("abs_min_pct","abs_min_time","abs_max_pct","abs_max_time",
            "rapid_min_pct","rapid_min_time","rapid_max_pct","rapid_max_time",
            "mid_mean_pct")
  missing <- setdiff(need, names(per_cell))
  if (length(missing) > 0) {
    stop("per_cell is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!("rapid_zero_inh" %in% names(per_cell))) per_cell$rapid_zero_inh <- FALSE
  if (!("rapid_zero_run_s" %in% names(per_cell))) per_cell$rapid_zero_run_s <- NA_real_


  # we expect these after the extract_bucket_features_one patch above
  if (!("mid_mean_pct_raw" %in% names(per_cell))) per_cell$mid_mean_pct_raw <- NA_real_
  if (!("mid_max_pct_raw"  %in% names(per_cell))) per_cell$mid_max_pct_raw  <- NA_real_

  has_area        <- "rapid_inh_area" %in% names(per_cell)
  has_bin_post_s  <- "bin_post_s" %in% names(per_cell)
  has_abs_persist <- has_cols("abs_has_inh_persist", "abs_has_exc_persist")

  vtrue <- function(x) dplyr::coalesce(as.logical(x), FALSE)

  per_cell %>%
    dplyr::mutate(
      has_rapid_exc = if (has_cols("rapid_has_exc_persist")) {
        vtrue(rapid_has_exc_persist)
      } else {
        is.finite(rapid_max_pct) & rapid_max_pct >= (100 + min_excursion_pct)
      },

      has_rapid_inh = if (has_cols("rapid_has_inh_persist")) {
        vtrue(rapid_has_inh_persist)
      } else {
        is.finite(rapid_min_pct) & rapid_min_pct <= (100 - min_excursion_pct)
      },

      is_biphasic_rapid = has_rapid_exc & has_rapid_inh,

      abs_exc_real = if (has_cols("max_is_real")) vtrue(max_is_real) else
        (is.finite(abs_max_pct) & abs_max_pct >= (100 + min_excursion_pct)),

      abs_inh_real = if (has_cols("min_is_real")) vtrue(min_is_real) else
        (is.finite(abs_min_pct) & abs_min_pct <= (100 - min_excursion_pct)),

      exc_mag = dplyr::if_else(is.finite(abs_max_pct), abs_max_pct - 100, NA_real_),
      inh_mag = dplyr::if_else(is.finite(abs_min_pct), 100 - abs_min_pct, NA_real_),

      inh_big_enough_rel =
        is.finite(inh_mag) & (
          (is.finite(exc_mag) & (inh_mag >= rel_inh_to_exc * exc_mag)) |
            (abs_inh_real & is.finite(abs_min_pct) & (abs_min_pct <= deep_dip_pct))
        ),

      inh_ok_time =
        is.finite(abs_min_time) &
        (abs_min_time >= rescue_inh_window[1]) &
        (abs_min_time <= rescue_inh_window[2]) &
        (if (has_bin_post_s) abs_min_time > pmax(0.1, bin_post_s) else abs_min_time > 0.1),

      rescue_inh_early =
        rescue_biphasic &
        inh_ok_time &
        is.finite(abs_min_pct) &
        (abs_min_pct <= (100 - (min_excursion_pct + rescue_inh_extra))),

      rescue_exc_later =
        rescue_biphasic &
        is.finite(abs_max_pct) & (abs_max_pct >= (100 + min_excursion_pct)) &
        is.finite(abs_max_time) & is.finite(abs_min_time) &
        (abs_max_time > abs_min_time + min_sep_s),

      dip_has_area = if (has_area) {
        is.finite(rapid_inh_area) & (rapid_inh_area >= rescue_area_thr)
      } else TRUE,

      dip_has_persist = if (has_abs_persist) {
        vtrue(abs_has_inh_persist) & vtrue(abs_has_exc_persist)
      } else TRUE,

      # ---- NEW: semi-smooth rebound support ----
      raw_mid_val = dplyr::case_when(
        raw_use == "mean" ~ mid_mean_pct_raw,
        raw_use == "max"  ~ mid_max_pct_raw,
        TRUE ~ mid_mean_pct_raw
      ),

      rebound_supported =
        is.finite(mid_mean_pct) & (mid_mean_pct >= (100 + rebound_mid_thr)) &
        is.finite(raw_mid_val)  & (raw_mid_val  >= (100 + raw_rebound_mid_thr)),

      not_still_inhibited_mid =
        is.finite(mid_mean_pct) & (mid_mean_pct >= (100 - still_inh_mid_thr)),

      is_biphasic_rescue =
        rescue_biphasic &
        rescue_inh_early & rescue_exc_later &
        dip_has_area & dip_has_persist &
        inh_big_enough_rel &
        not_still_inhibited_mid &
        rebound_supported,

      is_biphasic_rapid = is_biphasic_rapid & inh_big_enough_rel,

      biphasic_order = dplyr::case_when(
        vtrue(rapid_zero_inh) & abs_exc_real ~ "inh_then_exc",
        !(is_biphasic_rapid | is_biphasic_rescue) ~ NA_character_,
        is_biphasic_rapid & is.finite(rapid_max_time) & is.finite(rapid_min_time) ~
          dplyr::if_else(rapid_max_time < rapid_min_time, "exc_then_inh", "inh_then_exc"),
        is_biphasic_rescue & is.finite(abs_max_time) & is.finite(abs_min_time) ~
          dplyr::if_else(abs_max_time < abs_min_time, "exc_then_inh", "inh_then_exc"),
        TRUE ~ "ambiguous"
      ),

      biphasic_first = dplyr::case_when(
        biphasic_order == "exc_then_inh" ~ "excitation_first",
        biphasic_order == "inh_then_exc" ~ "inhibition_first",
        TRUE ~ NA_character_
      ),

      response_type = dplyr::case_when(
        # OVERRIDE: silence in rapid window is strong inhibition
        vtrue(rapid_zero_inh) & abs_exc_real ~ "biphasic",
        vtrue(rapid_zero_inh) & !abs_exc_real ~ "inhibition",

        # existing rules
        is_biphasic_rapid  ~ "biphasic",
        is_biphasic_rescue ~ "biphasic",
        abs_exc_real & abs_inh_real ~ "biphasic",
        abs_exc_real & !abs_inh_real ~ "excitation",
        abs_inh_real & !abs_exc_real ~ "inhibition",
        TRUE ~ "no_change"
      ),

    ) %>%
    dplyr::select(-dplyr::any_of("raw_mid_val"))
}

add_rapid_zero_inh <- function(per_cell,
                               df_prepped,
                               rapid_window = c(0, 10),
                               rate_col = "instRate",
                               zero_thr_hz = 0.1,
                               min_zero_s = 3) {

  stopifnot(all(c("cell_name", "time", rate_col) %in% names(df_prepped)))

  longest_true_run_seconds <- function(time, flag) {
    ok <- is.finite(time) & !is.na(flag)
    time <- time[ok]
    flag <- flag[ok]
    if (length(time) < 2) return(0)

    o <- order(time)
    time <- time[o]
    flag <- flag[o]

    # robust dt (for 1-sample runs)
    dt <- stats::median(diff(time), na.rm = TRUE)
    if (!is.finite(dt) || dt <= 0) dt <- 0

    r <- rle(flag)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    true_runs <- which(r$values)
    if (length(true_runs) == 0) return(0)

    durs <- vapply(true_runs, function(k) {
      i0 <- starts[k]
      i1 <- ends[k]
      (time[i1] - time[i0]) + dt
    }, numeric(1))

    max(durs, na.rm = TRUE)
  }

  zero_tbl <- df_prepped %>%
    dplyr::transmute(
      cell_name = .data$cell_name,
      time      = as.numeric(.data$time),
      rate      = as.numeric(.data[[rate_col]])
    ) %>%
    dplyr::filter(.data$time >= rapid_window[1], .data$time <= rapid_window[2]) %>%
    dplyr::filter(is.finite(.data$time), is.finite(.data$rate)) %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::summarise(
      rapid_zero_dur_s = longest_true_run_seconds(time, rate <= zero_thr_hz),
      rapid_zero_inh   = is.finite(rapid_zero_dur_s) & (rapid_zero_dur_s >= min_zero_s),
      .groups = "drop"
    )

  per_cell %>%
    dplyr::left_join(zero_tbl, by = "cell_name") %>%
    dplyr::mutate(
      rapid_zero_dur_s = dplyr::coalesce(.data$rapid_zero_dur_s, 0),
      rapid_zero_inh   = dplyr::coalesce(.data$rapid_zero_inh, FALSE)
    )
}


add_biphasic_features_simple <- function(
    per_cell,
    min_excursion_pct = 15,
    min_sep_s = 0.25,
    use_mid_rebound_gate = TRUE,
    rebound_mid_thr = 5,
    still_inh_mid_thr = 10
) {
  need <- c("abs_min_pct","abs_min_time","abs_max_pct","abs_max_time",
            "rapid_min_pct","rapid_min_time","rapid_max_pct","rapid_max_time")
  missing <- setdiff(need, names(per_cell))
  if (length(missing) > 0) {
    stop("per_cell is missing required columns: ", paste(missing, collapse = ", "))
  }

  # helper for safe logical flags
  vtrue <- function(x) dplyr::coalesce(as.logical(x), FALSE)

  per_cell %>%
    dplyr::mutate(
      # Presence of inh/exc evidence (prefer persist flags if you have them)
      has_rapid_inh = if ("rapid_has_inh_persist" %in% names(.)) vtrue(rapid_has_inh_persist) else {
        is.finite(rapid_min_pct) & (rapid_min_pct <= (100 - min_excursion_pct))
      },
      has_rapid_exc = if ("rapid_has_exc_persist" %in% names(.)) vtrue(rapid_has_exc_persist) else {
        is.finite(rapid_max_pct) & (rapid_max_pct >= (100 + min_excursion_pct))
      },

      # If you computed rapid_zero_inh, it should override as "strong inh"
      has_rapid_inh = has_rapid_inh | (if ("rapid_zero_inh" %in% names(.)) vtrue(rapid_zero_inh) else FALSE),

      # Absolute "real" (use your real flags if present)
      abs_inh_real = if ("min_is_real" %in% names(.)) vtrue(min_is_real) else {
        is.finite(abs_min_pct) & (abs_min_pct <= (100 - min_excursion_pct))
      },
      abs_exc_real = if ("max_is_real" %in% names(.)) vtrue(max_is_real) else {
        is.finite(abs_max_pct) & (abs_max_pct >= (100 + min_excursion_pct))
      },

      # Biphasic candidate (use absolute times to define order; enforce min_sep_s)
      is_biphasic = abs_inh_real & abs_exc_real &
        is.finite(abs_min_time) & is.finite(abs_max_time) &
        (abs(abs_max_time - abs_min_time) >= min_sep_s),

      biphasic_order = dplyr::case_when(
        !is_biphasic ~ NA_character_,
        abs_max_time < abs_min_time ~ "exc_then_inh",
        abs_min_time < abs_max_time ~ "inh_then_exc",
        TRUE ~ "ambiguous"
      ),

      # Optional mid gate: require rebound above baseline and not still inhibited
      mid_ok = if (use_mid_rebound_gate && ("mid_mean_pct" %in% names(.))) {
        is.finite(mid_mean_pct) &
          (mid_mean_pct >= (100 + rebound_mid_thr)) &
          (mid_mean_pct >= (100 - still_inh_mid_thr))
      } else {
        TRUE
      },

      is_biphasic = is_biphasic & mid_ok,

      # Final response_type
      response_type = dplyr::case_when(
        is_biphasic ~ "biphasic",
        abs_exc_real & !abs_inh_real ~ "excitation",
        abs_inh_real & !abs_exc_real ~ "inhibition",
        TRUE ~ "no_change"
      )
    )
}

#
#
# add_biphasic_features_simple <- function(
#     per_cell,
#     min_excursion_pct = 15,
#     min_sep_s = 0.25,
#     use_mid_rebound_gate = TRUE,
#     rebound_mid_thr = 5,
#     still_inh_mid_thr = 10
# ) {
#   has_cols <- function(...) all(c(...) %in% names(per_cell))
#   vtrue <- function(x) dplyr::coalesce(as.logical(x), FALSE)
#
#   # guard against accidental column name collision
#   if ("use_mid_rebound_gate" %in% names(per_cell)) {
#     stop("per_cell already has a column named `use_mid_rebound_gate` (name collision). Rename it or remove it.")
#   }
#
#   need <- c("abs_min_pct","abs_min_time","abs_max_pct","abs_max_time",
#             "rapid_min_pct","rapid_min_time","rapid_max_pct","rapid_max_time")
#   missing <- setdiff(need, names(per_cell))
#   if (length(missing) > 0) stop("per_cell missing: ", paste(missing, collapse = ", "))
#
#   # ensure expected optional columns exist
#   if (!("rapid_zero_inh" %in% names(per_cell))) per_cell$rapid_zero_inh <- FALSE
#   if (!("min_is_real" %in% names(per_cell)))     per_cell$min_is_real <- NA
#   if (!("max_is_real" %in% names(per_cell)))     per_cell$max_is_real <- NA
#   if (!("mid_mean_pct" %in% names(per_cell)))    per_cell$mid_mean_pct <- NA_real_
#
#   gate_mid <- isTRUE(use_mid_rebound_gate)  # force scalar TRUE/FALSE
#
#   out <- per_cell %>%
#     dplyr::mutate(
#       # --- canonical “real” flags ---
#       abs_inh_real = if (has_cols("min_is_real")) vtrue(min_is_real) else
#         (is.finite(abs_min_pct) & is.finite(abs_min_time) & abs_min_pct <= (100 - min_excursion_pct)),
#
#       abs_exc_real = if (has_cols("max_is_real")) vtrue(max_is_real) else
#         (is.finite(abs_max_pct) & is.finite(abs_max_time) & abs_max_pct >= (100 + min_excursion_pct)),
#
#       # silence in rapid window => strong inhibition
#       strong_rapid_inh = vtrue(rapid_zero_inh) | abs_inh_real,
#       strong_exc       = abs_exc_real,
#
#       # separation only matters if both exist
#       sep_ok = dplyr::if_else(
#         strong_rapid_inh & strong_exc & is.finite(abs_min_time) & is.finite(abs_max_time),
#         abs(abs_max_time - abs_min_time) >= min_sep_s,
#         TRUE
#       ),
#
#       # base classification
#       response_type = dplyr::case_when(
#         strong_rapid_inh & strong_exc & sep_ok ~ "biphasic",
#         strong_rapid_inh & !strong_exc ~ "inhibition",
#         !strong_rapid_inh & strong_exc ~ "excitation",
#         TRUE ~ "no_change"
#       ),
#
#       # order
#       biphasic_order = dplyr::case_when(
#         response_type != "biphasic" ~ NA_character_,
#         is.finite(abs_min_time) & is.finite(abs_max_time) ~
#           dplyr::if_else(abs_max_time > abs_min_time, "inh_then_exc", "exc_then_inh"),
#         is.finite(rapid_min_time) & is.finite(rapid_max_time) ~
#           dplyr::if_else(rapid_max_time > rapid_min_time, "inh_then_exc", "exc_then_inh"),
#         TRUE ~ "ambiguous"
#       ),
#
#       biphasic_first = dplyr::case_when(
#         response_type != "biphasic" ~ NA_character_,
#         biphasic_order == "exc_then_inh" ~ "excitation_first",
#         biphasic_order == "inh_then_exc" ~ "inhibition_first",
#         TRUE ~ NA_character_
#       )
#     )
#
#   # --- apply mid rebound gate OUTSIDE mutate (scalar if/else) ---
#   if (gate_mid) {
#     out <- out %>%
#       dplyr::mutate(
#         mid_ok     = is.finite(mid_mean_pct) & (mid_mean_pct >= (100 - still_inh_mid_thr)),
#         rebound_ok = is.finite(mid_mean_pct) & (mid_mean_pct >= (100 + rebound_mid_thr)),
#         response_type = dplyr::case_when(
#           response_type == "biphasic" &
#             biphasic_order == "inh_then_exc" &
#             (!mid_ok | !rebound_ok) ~ "inhibition",
#           TRUE ~ response_type
#         ),
#         # if we demoted it, clear order/first to avoid downstream “biphasic NA” weirdness
#         biphasic_order = dplyr::if_else(response_type == "biphasic", biphasic_order, NA_character_),
#         biphasic_first = dplyr::if_else(response_type == "biphasic", biphasic_first, NA_character_)
#       )
#   }
#
#   out
# }
#
#
#
# # Example: within each subclass, among biphasic cells only,
# # show exc-first vs inh-first proportions
# plot_biphasic_first_by_subclass <- function(per_cell2) {
#   per_cell2 %>%
#     dplyr::filter(response_type == "biphasic", !is.na(biphasic_first)) %>%
#     ggplot(aes(x = assigned_subclass, fill = biphasic_first)) +
#     geom_bar(position = "fill", color = "white") +
#     scale_y_continuous(labels = scales::percent_format()) +
#     labs(x = NULL, y = "Proportion of biphasic cells",
#          fill = "Biphasic order",
#          title = "Rapid biphasic order by subclass") +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank())
# }
# # show exc-first vs inh-first proportions
# plot_basic_by_subclass <- function(per_cell2) {
#   per_cell2 %>%
#     dplyr::filter(response_type != "biphasic") %>%
#     ggplot(aes(x = assigned_subclass, fill = response_type)) +
#     geom_bar(position = "fill", color = "white") +
#     scale_y_continuous(labels = scales::percent_format()) +
#     labs(x = NULL, y = "Proportion of biphasic cells",
#          fill = "Biphasic order",
#          title = "Rapid biphasic order by subclass") +
#     theme_bw() +
#     theme(panel.grid.minor = element_blank())
# }
#
