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
## 2) Persistence helper
## -----------------------------
# returns TRUE if there exists a continuous run (in time) beyond threshold
# of at least min_persist_s seconds inside the window.
has_persistent_excursion <- function(d, window, direction = c("exc","inh"),
                                     thr = 10,
                                     col = "pct_smooth",
                                     min_persist_s = 0.5) {
  direction <- match.arg(direction)

  d0 <- d %>%
    dplyr::filter(time >= window[1], time <= window[2]) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(time))

  if (nrow(d0) < 2) return(FALSE)

  flag <- if (direction == "exc") d0[[col]] >= (100 + thr) else d0[[col]] <= (100 - thr)

  if (!any(flag)) return(FALSE)

  # run-length encode contiguous TRUE segments (by index)
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  true_runs <- which(r$values)
  if (length(true_runs) == 0) return(FALSE)

  # compute duration of each TRUE run using time at segment boundaries
  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]
    dur <- d0$time[i1] - d0$time[i0]
    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }

  FALSE
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
                                        # NEW:
                                        min_persist_s = 0.5) {

  smooth_method <- match.arg(smooth_method)

  df_cell <- df_cell %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::arrange(time)

  df_cell <- add_smoothed_trace(df_cell,
                                k_smooth = k_smooth,
                                smooth_method = smooth_method,
                                loess_span = loess_span)

  get_extreme <- function(d, window, which = c("min","max"), col = "pct_smooth") {
    which <- match.arg(which)
    d0 <- d %>% dplyr::filter(time >= window[1], time <= window[2])
    v <- d0[[col]]
    if (nrow(d0) == 0 || all(!is.finite(v))) return(list(time = NA_real_, value = NA_real_))
    idx <- if (which == "min") which.min(v) else which.max(v)
    list(time = d0$time[idx], value = v[idx])
  }

  # --- use SMOOTH for decision-making extrema ---
  rapid_min <- get_extreme(df_cell, rapid_window, "min", col = "pct_smooth")
  rapid_max <- get_extreme(df_cell, rapid_window, "max", col = "pct_smooth")
  abs_min   <- get_extreme(df_cell, post_window,  "min", col = "pct_smooth")
  abs_max   <- get_extreme(df_cell, post_window,  "max", col = "pct_smooth")

  # --- persistence gates (smooth must stay beyond threshold for >= min_persist_s) ---
  rapid_has_exc_persist <- has_persistent_excursion(df_cell, rapid_window, "exc",
                                                    thr = min_excursion_pct,
                                                    col = "pct_smooth",
                                                    min_persist_s = min_persist_s)
  rapid_has_inh_persist <- has_persistent_excursion(df_cell, rapid_window, "inh",
                                                    thr = min_excursion_pct,
                                                    col = "pct_smooth",
                                                    min_persist_s = min_persist_s)

  # if you also want "absolute" excursions to require persistence, do it here:
  abs_has_exc_persist <- has_persistent_excursion(df_cell, post_window, "exc",
                                                  thr = min_excursion_pct,
                                                  col = "pct_smooth",
                                                  min_persist_s = min_persist_s)
  abs_has_inh_persist <- has_persistent_excursion(df_cell, post_window, "inh",
                                                  thr = min_excursion_pct,
                                                  col = "pct_smooth",
                                                  min_persist_s = min_persist_s)

  # Flag weak extrema using BOTH threshold + persistence
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

  # Rapid biphasic logic uses PERSISTENT detections (not single points)
  rapid_pattern <- dplyr::case_when(
    rapid_has_exc_persist & rapid_has_inh_persist ~
      ifelse(rapid_max$time < rapid_min$time, "rapid_biphasic_exc_to_inh", "rapid_biphasic_inh_to_exc"),
    rapid_has_exc_persist ~ "rapid_excitation_only",
    rapid_has_inh_persist ~ "rapid_inhibition_only",
    TRUE ~ "rapid_no_clear_change"
  )

  # Mid/late summaries still from smooth
  mid_mean <- df_cell %>%
    dplyr::filter(time >= mid_window[1], time <= mid_window[2]) %>%
    dplyr::summarise(x = mean(pct_smooth, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(x)

  late_mean <- df_cell %>%
    dplyr::filter(time >= late_mean_window[1], time <= late_mean_window[2]) %>%
    dplyr::summarise(x = mean(pct_smooth, na.rm = TRUE), .groups = "drop") %>%
    dplyr::pull(x)

  late_slope <- {
    dlate <- df_cell %>%
      dplyr::filter(time >= late_slope_window[1], time <= late_slope_window[2]) %>%
      dplyr::filter(is.finite(time), is.finite(pct_smooth))
    if (nrow(dlate) >= 5) coef(stats::lm(pct_smooth ~ time, data = dlate))[["time"]] else NA_real_
  }

  tibble::tibble(
    cell_name = df_cell$cell_name[1],
    assigned_subclass = if ("assigned_subclass" %in% names(df_cell)) df_cell$assigned_subclass[1] else NA_character_,

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



extract_bucket_features <- function(df_prepped,
                                    rapid_window = c(0, 10),
                                    mid_window   = c(10, 40),
                                    late_slope_window = c(40, 60),
                                    late_mean_window  = c(50, 60),
                                    post_window = c(0, 60),
                                    min_excursion_pct = 5,
                                    k_smooth = 11,
                                    smooth_method = c("rollmedian","rollmean","loess"),
                                    loess_span = 0.15) {

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
      loess_span = loess_span
    )) %>%
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
    { if (isTRUE(feats$max_is_real) && is.finite(feats$abs_max_time) && is.finite(feats$abs_max_pct_raw))
      geom_point(data = tibble::tibble(time = feats$abs_max_time, y = feats$abs_max_pct_raw),
                 aes(x = time, y = y), color = abs_max_col, size = 2.6)
      else NULL } +
    { if (isTRUE(feats$min_is_real) && is.finite(feats$abs_min_time) && is.finite(feats$abs_min_pct_raw))
      geom_point(data = tibble::tibble(time = feats$abs_min_time, y = feats$abs_min_pct_raw),
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
add_biphasic_features <- function(per_cell, min_excursion_pct = 15) {

  need <- c("rapid_min_pct","rapid_max_pct","rapid_min_time","rapid_max_time",
            "abs_min_pct","abs_max_pct")
  missing <- setdiff(need, names(per_cell))
  if (length(missing) > 0) {
    stop("per_cell is missing required columns: ", paste(missing, collapse = ", "))
  }

  per_cell %>%
    dplyr::mutate(
      # thresholds relative to baseline (=100)
      has_rapid_exc = is.finite(rapid_max_pct) & rapid_max_pct >= (100 + min_excursion_pct),
      has_rapid_inh = is.finite(rapid_min_pct) & rapid_min_pct <= (100 - min_excursion_pct),

      is_biphasic_rapid = has_rapid_exc & has_rapid_inh,

      # order of rapid biphasic extrema
      biphasic_order = dplyr::case_when(
        !is_biphasic_rapid ~ NA_character_,
        !is.finite(rapid_max_time) | !is.finite(rapid_min_time) ~ "ambiguous",
        rapid_max_time < rapid_min_time ~ "exc_then_inh",
        rapid_min_time < rapid_max_time ~ "inh_then_exc",
        TRUE ~ "ambiguous"
      ),

      # optional convenience label: what came first
      biphasic_first = dplyr::case_when(
        biphasic_order == "exc_then_inh" ~ "excitation_first",
        biphasic_order == "inh_then_exc" ~ "inhibition_first",
        TRUE ~ NA_character_
      ),

      # response type (biphasic takes precedence)
      response_type = dplyr::case_when(
        is_biphasic_rapid ~ "biphasic",
        is.finite(abs_max_pct) & abs_max_pct >= (100 + min_excursion_pct) ~ "excitation",
        is.finite(abs_min_pct) & abs_min_pct <= (100 - min_excursion_pct) ~ "inhibition",
        TRUE ~ "no_change"
      )
    )
}

add_biphasic_features_consistent <- function(per_cell, min_excursion_pct = 15) {

  # if you have these columns, we will use them; otherwise we fall back safely
  has_cols <- function(...) all(c(...) %in% names(per_cell))

  per_cell %>%
    dplyr::mutate(
      # Prefer persistence-gated truth if present
      has_rapid_exc = if (has_cols("rapid_has_exc_persist")) rapid_has_exc_persist
      else is.finite(rapid_max_pct) & rapid_max_pct >= (100 + min_excursion_pct),

      has_rapid_inh = if (has_cols("rapid_has_inh_persist")) rapid_has_inh_persist
      else is.finite(rapid_min_pct) & rapid_min_pct <= (100 - min_excursion_pct),

      is_biphasic_rapid = has_rapid_exc & has_rapid_inh,

      biphasic_order = dplyr::case_when(
        !is_biphasic_rapid ~ NA_character_,
        !is.finite(rapid_max_time) | !is.finite(rapid_min_time) ~ "ambiguous",
        rapid_max_time < rapid_min_time ~ "exc_then_inh",
        rapid_min_time < rapid_max_time ~ "inh_then_exc",
        TRUE ~ "ambiguous"
      ),

      biphasic_first = dplyr::case_when(
        biphasic_order == "exc_then_inh" ~ "excitation_first",
        biphasic_order == "inh_then_exc" ~ "inhibition_first",
        TRUE ~ NA_character_
      ),

      # Absolute “real” flags should drive mono-phasic calls
      abs_exc_real = if (has_cols("max_is_real")) max_is_real
      else is.finite(abs_max_pct) & abs_max_pct >= (100 + min_excursion_pct),

      abs_inh_real = if (has_cols("min_is_real")) min_is_real
      else is.finite(abs_min_pct) & abs_min_pct <= (100 - min_excursion_pct),

      # response_type: biphasic takes precedence
      response_type = dplyr::case_when(
        is_biphasic_rapid ~ "biphasic",
        abs_exc_real & !abs_inh_real ~ "excitation",
        abs_inh_real & !abs_exc_real ~ "inhibition",
        abs_exc_real & abs_inh_real ~ "biphasic",   # optional: classify abs-biphasic as biphasic
        TRUE ~ "no_change"
      )
    )
}

# Example: within each subclass, among biphasic cells only,
# show exc-first vs inh-first proportions
plot_biphasic_first_by_subclass <- function(per_cell2) {
  per_cell2 %>%
    dplyr::filter(response_type == "biphasic", !is.na(biphasic_first)) %>%
    ggplot(aes(x = assigned_subclass, fill = biphasic_first)) +
    geom_bar(position = "fill", color = "white") +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = NULL, y = "Proportion of biphasic cells",
         fill = "Biphasic order",
         title = "Rapid biphasic order by subclass") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
}

# =========================================================
# Your requested factor re-leveling (UPDATED)
# =========================================================
## Adjust here for excursion_pct boundaries
# per_cell2 <- per_cell %>%
#   add_biphasic_features(min_excursion_pct = 15) %>%
#   dplyr::mutate(
#     biphasic_order = factor(
#       biphasic_order,
#       levels = c("exc_then_inh","inh_then_exc","ambiguous")
#     ),
#     response_type = factor(
#       response_type,
#       levels = c("excitation","inhibition","biphasic","no_change")
#     )
#   )
