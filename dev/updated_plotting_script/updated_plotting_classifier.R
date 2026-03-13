# puff_qc.R
# Consolidated functions for puff-response QC, feature extraction, and plotting
# Author: (refactor for Scott S.)
# Usage: source("puff_qc.R"); then call prep_puff_df(), extract_bucket_features(), plot_bucket_qc(), etc.

# ---- packages ----------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(zoo)
  library(tibble)
  library(scales)
  library(forcats)
  library(stringr)
})

# ---- 0. small helpers ---------------------------------------------------------

# simple centered moving average (symmetric filter)
smooth_ma <- function(x, k = 7) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

# linear fit summary (slope, r2, p)
linfit_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(slope = NA_real_, r2 = NA_real_, p = NA_real_))
  fit <- stats::lm(y ~ x)
  s <- summary(fit)
  list(
    slope = unname(coef(fit)[2]),
    r2    = s$r.squared,
    p     = s$coefficients[2, 4]
  )
}

# run-length contiguous-duration test for "persistence" of excursion
has_persistent_excursion <- function(d, window, direction = c("exc","inh"),
                                     thr = 10, col = "pct_smooth",
                                     min_persist_s = 0.5) {
  direction <- match.arg(direction)
  d0 <- d %>%
    dplyr::filter(.data$time >= window[1], .data$time <= window[2]) %>%
    dplyr::arrange(.data$time) %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(.data$time))
  if (nrow(d0) < 2) return(FALSE)
  flag <- if (direction == "exc") d0[[col]] >= (100 + thr) else d0[[col]] <= (100 - thr)
  if (!any(flag)) return(FALSE)
  r <- rle(flag)
  ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  true_runs <- which(r$values)
  if (length(true_runs) == 0) return(FALSE)
  for (k in true_runs) {
    i0 <- starts[k]; i1 <- ends[k]
    dur <- d0$time[i1] - d0$time[i0]
    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }
  FALSE
}

# ---- 1. preprocessing ---------------------------------------------------------

#' Prepare long df for puff analysis
#'
#' Expects columns: cell_name, time, instRate (others allowed).
#' Adds per-cell baseline_hz and pct_of_baseline (100 = baseline).
#'
#' @param df data.frame long format (timepoints)
#' @param baseline_window numeric length-2 baseline window (seconds)
#' @return data.frame with baseline_hz, pct_of_baseline
prep_puff_df <- function(df, baseline_window = c(-5, 0)) {
  stopifnot(is.data.frame(df))
  stopifnot(all(c("cell_name", "time", "instRate") %in% names(df)))

  df0 <- df %>%
    dplyr::mutate(
      time = as.numeric(.data$time),
      instRate = as.numeric(.data$instRate)
    )

  base_tbl <- df0 %>%
    dplyr::filter(.data$time >= baseline_window[1], .data$time < baseline_window[2]) %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::summarise(
      baseline_hz = mean(.data$instRate, na.rm = TRUE),
      n_baseline_pts = sum(is.finite(.data$instRate)),
      .groups = "drop"
    )

  df1 <- df0 %>%
    dplyr::left_join(base_tbl, by = "cell_name") %>%
    dplyr::mutate(
      pct_of_baseline = ifelse(is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
                               100 * .data$instRate / .data$baseline_hz,
                               NA_real_)
    )

  df1
}

# ---- 2. smoothing utilities ---------------------------------------------------

#' Add smoothed pct-of-baseline trace to a single-cell df
#' supports rollmedian, rollmean, or loess smoothing
add_smoothed_trace <- function(df_cell,
                               k_smooth = 11,
                               smooth_method = c("rollmedian", "rollmean", "loess"),
                               loess_span = 0.15) {
  smooth_method <- match.arg(smooth_method)
  df_cell <- df_cell %>% dplyr::arrange(.data$time) %>% dplyr::mutate(time = as.numeric(.data$time))
  y <- df_cell$pct_of_baseline

  y_smooth <- switch(smooth_method,
                     rollmedian = zoo::rollmedian(y, k = k_smooth, fill = NA, align = "center"),
                     rollmean   = zoo::rollmean(y,   k = k_smooth, fill = NA, align = "center"),
                     loess = {
                       ok <- is.finite(df_cell$time) & is.finite(y)
                       if (sum(ok) < 10) rep(NA_real_, length(y)) else {
                         fit <- stats::loess(y[ok] ~ df_cell$time[ok], span = loess_span, degree = 1)
                         out <- rep(NA_real_, length(y)); out[ok] <- stats::predict(fit, newdata = df_cell$time[ok]); out
                       }
                     }
  )

  df_cell %>% dplyr::mutate(pct_smooth = as.numeric(y_smooth))
}

# ---- 3. per-cell bucketed feature extraction ---------------------------------

#' Extract bucketed features from a single cell long df
#'
#' Returns a 1-row tibble of features: rapid/min/max timings, persistence flags,
#' mid/late means/slopes, direction, peak magnitudes, and meta about smoothing.
#'
#' @param df_cell long df (single cell) with pct_of_baseline
#' @param rapid_window, mid_window, late_slope_window, late_mean_window, post_window numeric windows
#' @param min_excursion_pct numeric threshold in percent from baseline (100)
#' @param k_smooth smoothing kernel size
#' @param smooth_method smoothing method ("rollmedian","rollmean","loess")
#' @param loess_span only for loess
#' @param min_persist_s minimal continuous seconds to call an excursion persistent
#' @return tibble 1 row (features)
extract_bucket_features_one <- function(df_cell,
                                        rapid_window = c(0, 10),
                                        mid_window   = c(10, 40),
                                        late_slope_window = c(40, 60),
                                        late_mean_window  = c(50, 60),
                                        post_window = c(0, 60),
                                        min_excursion_pct = 5,
                                        k_smooth = 11,
                                        smooth_method = c("rollmedian", "rollmean", "loess"),
                                        loess_span = 0.15,
                                        min_persist_s = 0.5) {

  smooth_method <- match.arg(smooth_method)
  df_cell <- df_cell %>% dplyr::mutate(time = as.numeric(.data$time)) %>% dplyr::arrange(.data$time)

  # add smoothed trace
  df_cell <- add_smoothed_trace(df_cell, k_smooth = k_smooth, smooth_method = smooth_method, loess_span = loess_span)

  get_extreme <- function(d, window, which = c("min", "max"), col = "pct_smooth") {
    which <- match.arg(which)
    d0 <- d %>% dplyr::filter(.data$time >= window[1], .data$time <= window[2])
    v <- d0[[col]]
    if (nrow(d0) == 0 || all(!is.finite(v))) return(list(time = NA_real_, value = NA_real_))
    idx <- if (which == "min") which.min(v) else which.max(v)
    list(time = d0$time[idx], value = v[idx])
  }

  rapid_min <- get_extreme(df_cell, rapid_window, "min", col = "pct_smooth")
  rapid_max <- get_extreme(df_cell, rapid_window, "max", col = "pct_smooth")
  abs_min   <- get_extreme(df_cell, post_window,  "min", col = "pct_smooth")
  abs_max   <- get_extreme(df_cell, post_window,  "max", col = "pct_smooth")

  rapid_has_exc_persist <- has_persistent_excursion(df_cell, rapid_window, "exc", thr = min_excursion_pct, col = "pct_smooth", min_persist_s = min_persist_s)
  rapid_has_inh_persist <- has_persistent_excursion(df_cell, rapid_window, "inh", thr = min_excursion_pct, col = "pct_smooth", min_persist_s = min_persist_s)

  abs_has_exc_persist <- has_persistent_excursion(df_cell, post_window, "exc", thr = min_excursion_pct, col = "pct_smooth", min_persist_s = min_persist_s)
  abs_has_inh_persist <- has_persistent_excursion(df_cell, post_window, "inh", thr = min_excursion_pct, col = "pct_smooth", min_persist_s = min_persist_s)

  max_is_real <- is.finite(abs_max$value) && abs_max$value >= (100 + min_excursion_pct) && abs_has_exc_persist
  min_is_real <- is.finite(abs_min$value) && abs_min$value <= (100 - min_excursion_pct) && abs_has_inh_persist

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
    rapid_has_exc_persist & rapid_has_inh_persist ~ ifelse(rapid_max$time < rapid_min$time, "rapid_biphasic_exc_to_inh", "rapid_biphasic_inh_to_exc"),
    rapid_has_exc_persist ~ "rapid_excitation_only",
    rapid_has_inh_persist ~ "rapid_inhibition_only",
    TRUE ~ "rapid_no_clear_change"
  )

  # mid mean and late slope
  mid_mean <- df_cell %>% dplyr::filter(.data$time >= mid_window[1], .data$time <= mid_window[2]) %>% dplyr::summarise(x = mean(.data$pct_smooth, na.rm = TRUE), .groups = "drop") %>% dplyr::pull(x)
  late_mean <- df_cell %>% dplyr::filter(.data$time >= late_mean_window[1], .data$time <= late_mean_window[2]) %>% dplyr::summarise(x = mean(.data$pct_smooth, na.rm = TRUE), .groups = "drop") %>% dplyr::pull(x)

  late_slope <- {
    dlate <- df_cell %>% dplyr::filter(.data$time >= late_slope_window[1], .data$time <= late_slope_window[2]) %>% dplyr::filter(is.finite(.data$time), is.finite(.data$pct_smooth))
    if (nrow(dlate) >= 5) coef(stats::lm(pct_smooth ~ time, data = dlate))[["time"]] else NA_real_
  }

  tibble::tibble(
    #cell_name = NA_character_,
    assigned_subclass = if ("assigned_subclass" %in% names(df_cell)) df_cell$assigned_subclass[1] else NA_character_,
    baseline_hz = if ("baseline_hz" %in% names(df_cell)) df_cell$baseline_hz[1] else NA_real_,

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

    k_smooth = k_smooth,
    smooth_method = smooth_method,
    min_persist_s = min_persist_s
  )
}

# wrapper to run across all cells
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
    dplyr::group_by(.data$cell_name) %>%
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

add_response_call <- function(per_cell,
                              mild = 10,
                              moderate = 25,
                              strong = 50) {

  stopifnot(is.data.frame(per_cell))
  need <- c("response_type", "exc_mag", "inh_mag")
  missing <- setdiff(need, names(per_cell))
  if (length(missing) > 0) stop("add_response_call: missing columns: ", paste(missing, collapse = ", "))

  classify_mag <- function(mag) {
    dplyr::case_when(
      !is.finite(mag) ~ NA_character_,
      mag >= strong ~ "strong",
      mag >= moderate ~ "moderate",
      mag >= mild ~ "mild",
      TRUE ~ "no_change"
    )
  }

  per_cell %>%
    dplyr::mutate(
      # pick the magnitude consistent with the call
      mag_for_call = dplyr::case_when(
        response_type == "excitation" ~ exc_mag,
        response_type == "inhibition" ~ inh_mag,
        TRUE ~ NA_real_
      ),
      mag_bin = classify_mag(mag_for_call),

      response_call = dplyr::case_when(
        response_type == "excitation" & mag_bin == "mild" ~ "mild_excitation",
        response_type == "excitation" & mag_bin == "moderate" ~ "moderate_excitation",
        response_type == "excitation" & mag_bin == "strong" ~ "strong_excitation",

        response_type == "inhibition" & mag_bin == "mild" ~ "mild_inhibition",
        response_type == "inhibition" & mag_bin == "moderate" ~ "moderate_inhibition",
        response_type == "inhibition" & mag_bin == "strong" ~ "strong_inhibition",

        response_type %in% c("no_change", "biphasic") ~ "no_change",
        TRUE ~ "no_change"
      ),
      response_call = factor(response_call, levels = response_levels)
    )
}


# ---- 4. QC single-cell plotting ----------------------------------------------

#' Plot QC for a single cell
#'
#' Requires df_prepped (from prep_puff_df) with pct_of_baseline.
plot_bucket_qc <- function(df_prepped,
                           cell_id,
                           rapid_window = c(0,10),
                           mid_window = c(10,40),
                           late_slope_window = c(40,60),
                           late_mean_window = c(50,60),
                           post_window = c(0,50),
                           k_smooth = 11,
                           smooth_method = c("rollmedian","rollmean","loess"),
                           loess_span = 0.15,
                           min_excursion_pct = 5,
                           show_instRate = FALSE) {

  smooth_method <- match.arg(smooth_method)

  d <- df_prepped %>%
    dplyr::filter(.data$cell_name == cell_id) %>%
    dplyr::mutate(time = as.numeric(.data$time)) %>%
    dplyr::arrange(.data$time)

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

  # late fit
  dlate <- d %>% dplyr::filter(.data$time >= late_slope_window[1], .data$time <= late_slope_window[2]) %>% dplyr::filter(is.finite(.data$time), is.finite(.data$pct_smooth))
  fit_line <- NULL
  if (nrow(dlate) >= 5) {
    fit <- stats::lm(pct_smooth ~ time, data = dlate)
    fit_line <- tibble::tibble(time = seq(min(dlate$time), max(dlate$time), length.out = 50)) %>% dplyr::mutate(pct_fit = stats::predict(fit, newdata = .))
  }

  col_exc <- "#d62728"  # excitation (red)
  col_inh <- "#1f77b4"  # inhibition (blue)

  title_txt <- paste0(cell_id,
                      if (!is.na(feats$assigned_subclass)) paste0(" | ", feats$assigned_subclass) else "",
                      if (!is.na(feats$direction)) paste0(" | ", feats$direction) else "",
                      if (!is.na(feats$response_peak)) paste0(" | peak = ", sprintf("%.1f", feats$response_peak), "%") else "")

  subtitle_txt <- paste0(
    "Baseline = ", ifelse(is.finite(feats$baseline_hz), sprintf("%.2f Hz", feats$baseline_hz), "NA"),
    " | rapid = ", feats$rapid_pattern,
    " | abs_min=", ifelse(is.finite(feats$abs_min_pct), sprintf("%.1f%%", feats$abs_min_pct), "NA/flagged"),
    " | abs_max=", ifelse(is.finite(feats$abs_max_pct), sprintf("%.1f%%", feats$abs_max_pct), "NA/flagged")
  )

  p <- ggplot(d, aes(x = time)) +
    annotate("rect", xmin = rapid_window[1], xmax = rapid_window[2], ymin = -Inf, ymax = Inf, alpha = 0.08) +
    annotate("rect", xmin = mid_window[1], xmax = mid_window[2], ymin = -Inf, ymax = Inf, alpha = 0.05) +
    annotate("rect", xmin = late_slope_window[1], xmax = late_slope_window[2], ymin = -Inf, ymax = Inf, alpha = 0.06) +
    geom_line(aes(y = pct_of_baseline), color = "grey60", linewidth = 0.3, alpha = 0.7) +
    geom_line(aes(y = pct_smooth), color = "black", linewidth = 0.7, na.rm = TRUE) +
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5, color = "grey30") +
    geom_hline(yintercept = 100 + min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_exc) +
    geom_hline(yintercept = 100 - min_excursion_pct, linetype = 3, linewidth = 0.4, color = col_inh) +
    { if (!is.null(fit_line)) geom_line(data = fit_line, aes(x = time, y = pct_fit), linewidth = 1.0) else NULL } +
    labs(title = title_txt, subtitle = subtitle_txt, x = "Time (s)", y = "% of baseline (100 = baseline)") +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

  p
}

# ---- 5. group plotting helpers / visualization --------------------------------

# balanced log transforms used in your originals (kept as utilities)
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

# assemble trace + mean data for EI plots
make_plot_df_ei <- function(df_prepped, per_cell2,
                            xlim = c(-5, 56), k_smooth = 3,
                            time_bin_s = NULL, min_cells_per_time = 3) {

  if (is.null(time_bin_s)) {
    dt_est <- df_prepped %>%
      dplyr::group_by(.data$cell_name) %>%
      dplyr::summarise(dt = stats::median(diff(sort(unique(.data$time)))), .groups = "drop") %>%
      dplyr::summarise(dt = stats::median(.data$dt, na.rm = TRUE)) %>%
      dplyr::pull(.data$dt)
    time_bin_s <- max(0.1, dt_est)
  }

  d0 <- df_prepped %>%
    dplyr::select(-dplyr::any_of("assigned_subclass")) %>%
    dplyr::inner_join(per_cell2 %>% dplyr::select(.data$cell_name, .data$response_type, .data$assigned_subclass), by = "cell_name") %>%
    dplyr::filter(.data$response_type %in% c("excitation", "inhibition"),
                  .data$time >= xlim[1], .data$time <= xlim[2],
                  is.finite(.data$pct_of_baseline)) %>%
    dplyr::mutate(time = as.numeric(.data$time))

  d_traces <- d0 %>%
    dplyr::group_by(.data$response_type, .data$assigned_subclass, .data$cell_name) %>%
    dplyr::arrange(.data$time, .by_group = TRUE) %>%
    dplyr::mutate(pct_smooth = zoo::rollmean(.data$pct_of_baseline, k = k_smooth, fill = NA, align = "center")) %>%
    dplyr::ungroup()

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

# ---- 6. summary / classification plots (stacked, heatmap, dots, donuts) -----

# helper: ensure assigned_subclass present (left-join from raw df if needed)
ensure_subclass <- function(per_cell, raw_df = NULL, cell_col = "cell_name") {
  stopifnot(is.data.frame(per_cell)); stopifnot(cell_col %in% names(per_cell))
  if ("assigned_subclass" %in% names(per_cell)) return(per_cell)
  if (is.null(raw_df)) stop("per_cell has no 'assigned_subclass'. Provide raw_df to join it.")
  stopifnot(is.data.frame(raw_df)); stopifnot(cell_col %in% names(raw_df)); stopifnot("assigned_subclass" %in% names(raw_df))
  map_sub <- raw_df %>% dplyr::select(dplyr::all_of(c(cell_col, "assigned_subclass"))) %>% dplyr::distinct()
  per_cell %>% dplyr::left_join(map_sub, by = setNames(cell_col, cell_col))
}

prep_subclasses <- function(per_cell, keep = c("L23_IT", "L5_ET", "L5_IT"), recode_map = c("L3c" = "L23_IT")) {
  stopifnot(is.data.frame(per_cell)); stopifnot("assigned_subclass" %in% names(per_cell))
  per_cell %>%
    dplyr::mutate(
      assigned_subclass = dplyr::if_else(.data$assigned_subclass %in% names(recode_map), unname(recode_map[.data$assigned_subclass]), .data$assigned_subclass),
      assigned_subclass = as.character(.data$assigned_subclass)
    ) %>%
    dplyr::filter(.data$assigned_subclass %in% keep) %>%
    dplyr::mutate(assigned_subclass = factor(.data$assigned_subclass, levels = keep))
}

# response color and label maps (kept from originals)
response_levels <- c("strong_inhibition","moderate_inhibition","mild_inhibition","no_change","mild_excitation","moderate_excitation","strong_excitation")
response_cols <- c(strong_inhibition = "#B2182B", moderate_inhibition = "#D6604D", mild_inhibition = "#F4A582",
                   no_change = "grey80", mild_excitation = "#92C5DE", moderate_excitation = "#4393C3", strong_excitation = "#2166AC")
pretty_labels <- setNames(c("Strong inhibition","Moderate inhibition","Mild inhibition","No change","Mild excitation","Moderate excitation","Strong excitation"), response_levels)

prep_plot_df <- function(per_cell) {
  stopifnot("response_call" %in% names(per_cell)); stopifnot("assigned_subclass" %in% names(per_cell))
  per_cell %>% dplyr::filter(!is.na(.data$response_call)) %>% dplyr::mutate(response_call = factor(.data$response_call, levels = response_levels))
}

plot_stacked_100 <- function(per_cell, title = "Puff response composition by subclass", show_counts_in_strip = TRUE) {
  dat <- prep_plot_df(per_cell)
  n_sub <- dat %>% dplyr::count(.data$assigned_subclass, name = "N")
  p <- dat %>% dplyr::count(.data$assigned_subclass, .data$response_call, name = "n") %>% dplyr::group_by(.data$assigned_subclass) %>% dplyr::mutate(frac = .data$n / sum(.data$n)) %>% dplyr::ungroup() %>%
    ggplot(aes(x = .data$assigned_subclass, y = .data$frac, fill = .data$response_call)) +
    geom_col(width = 0.72, color = "white", linewidth = 0.3) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_fill_manual(values = response_cols, drop = FALSE, labels = pretty_labels) +
    labs(x = NULL, y = "Percent of cells", fill = "Response class", title = title) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(), legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))
  if (show_counts_in_strip) p <- p + geom_text(data = n_sub, aes(x = .data$assigned_subclass, y = 1.02, label = paste0("N=", .data$N)), inherit.aes = FALSE, size = 3.6, vjust = 0) + coord_cartesian(clip = "off")
  p
}

plot_heatmap_prop <- function(per_cell, title = "Response proportions (heatmap)") {
  dat <- prep_plot_df(per_cell)
  hm <- dat %>% dplyr::count(.data$assigned_subclass, .data$response_call, name = "n") %>% dplyr::group_by(.data$assigned_subclass) %>% dplyr::mutate(frac = .data$n / sum(.data$n)) %>% dplyr::ungroup()
  ggplot(hm, aes(x = .data$response_call, y = .data$assigned_subclass, fill = .data$frac)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = percent(.data$frac, accuracy = 1)), size = 3.4) +
    scale_x_discrete(labels = pretty_labels) + scale_y_discrete(drop = FALSE) +
    scale_fill_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = NULL, y = NULL, fill = "Percent", title = title) +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 35, hjust = 1), plot.title = element_text(face = "bold", hjust = 0.5))
}

plot_prop_points <- function(per_cell, title = "Response proportions with 95% CI") {
  dat <- prep_plot_df(per_cell)
  dfp <- dat %>% dplyr::count(.data$assigned_subclass, .data$response_call, name = "n") %>% dplyr::group_by(.data$assigned_subclass) %>% dplyr::mutate(N = sum(.data$n), frac = .data$n / .data$N) %>% dplyr::ungroup()
  dfp <- dfp %>% dplyr::mutate(se = sqrt(pmax(frac * (1 - frac) / pmax(N, 1), 0)), lo = pmax(frac - 1.96 * se, 0), hi = pmin(frac + 1.96 * se, 1))
  ggplot(dfp, aes(x = .data$response_call, y = .data$frac, color = .data$response_call)) +
    geom_point(size = 2.6) +
    geom_errorbar(aes(ymin = .data$lo, ymax = .data$hi), width = 0.15, linewidth = 0.5) +
    facet_wrap(~assigned_subclass, nrow = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    scale_x_discrete(labels = pretty_labels) +
    scale_color_manual(values = response_cols, drop = FALSE, guide = "none") +
    labs(x = NULL, y = "Percent of cells", title = title) +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))
}

# small convenience: per-cell violin with baseline-centered y axis
plot_centered_violin <- function(per_cell, title = "Per-cell response magnitude (baseline-centered QC)") {
  stopifnot("assigned_subclass" %in% names(per_cell)); stopifnot("response_peak" %in% names(per_cell))
  ggplot(per_cell, aes(x = .data$assigned_subclass, y = .data$response_peak)) +
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.6, color = "black") +
    geom_violin(fill = "grey85", color = "black", linewidth = 0.4, trim = TRUE) +
    geom_point(aes(color = .data$direction), position = position_jitter(width = 0.12), size = 1.4, alpha = 0.6) +
    scale_color_manual(values = c(inhibition = "#B2182B", excitation = "#2166AC"), guide = "none") +
    coord_flip() + labs(x = NULL, y = "Percent of baseline (100 = no change)", title = title) +
    theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold", hjust = 0.5), panel.grid.minor = element_blank())
}

# ---- 7. batch helpers ---------------------------------------------------------

plot_bucket_qc_many <- function(df_prepped, cell_ids, out_dir = "qc_bucket_plots", width = 7, height = 4.5, dpi = 200, ...) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (cid in cell_ids) {
    p <- try(plot_bucket_qc(df_prepped, cid, ...), silent = TRUE)
    if (!inherits(p, "try-error")) {
      fn <- file.path(out_dir, paste0(gsub("[^A-Za-z0-9_.-]", "_", cid), ".png"))
      ggsave(fn, p, width = width, height = height, dpi = dpi)
    }
  }
  message("Saved ", length(cell_ids), " QC plots to: ", out_dir)
}


add_biphasic_features <- function(per_cell, min_excursion_pct = 5) {
  need <- c("rapid_has_exc_persist","rapid_has_inh_persist","min_is_real","max_is_real",
            "rapid_min_time","rapid_max_time")
  missing <- setdiff(need, names(per_cell))
  if (length(missing) > 0) stop("add_biphasic_features: missing: ", paste(missing, collapse = ", "))

  per_cell %>%
    dplyr::mutate(
      has_rapid_exc = .data$rapid_has_exc_persist,
      has_rapid_inh = .data$rapid_has_inh_persist,
      is_biphasic_rapid = .data$has_rapid_exc & .data$has_rapid_inh,

      biphasic_order = dplyr::case_when(
        !.data$is_biphasic_rapid ~ NA_character_,
        !is.finite(.data$rapid_max_time) | !is.finite(.data$rapid_min_time) ~ "ambiguous",
        .data$rapid_max_time < .data$rapid_min_time ~ "exc_then_inh",
        .data$rapid_min_time < .data$rapid_max_time ~ "inh_then_exc",
        TRUE ~ "ambiguous"
      ),

      biphasic_first = dplyr::case_when(
        biphasic_order == "exc_then_inh" ~ "excitation_first",
        biphasic_order == "inh_then_exc" ~ "inhibition_first",
        TRUE ~ NA_character_
      ),

      response_type = dplyr::case_when(
        .data$is_biphasic_rapid ~ "biphasic",
        .data$max_is_real & !.data$min_is_real ~ "excitation",
        .data$min_is_real & !.data$max_is_real ~ "inhibition",
        .data$min_is_real & .data$max_is_real ~ "biphasic",
        TRUE ~ "no_change"
      )
    )
}

# ---- 8. small examples / workflow snippet -----------------------------------
# Example workflow (uncomment & adapt in your script):


