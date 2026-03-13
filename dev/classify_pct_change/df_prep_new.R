## =========================================================
## Puff response QC + bucketed feature extraction + (prep) tables
##
## CLEAN BUILD SCRIPT (copy/paste)
##
## Assumes df has columns:
##   time, instRate, cell_name, assigned_subclass
## Optional but supported:
##   puff, bath, blockers, protocol
##
## Key rules:
## - Baseline = mean(instRate) in [-5, 0) seconds (per cell)
## - pct_of_baseline = 100 * instRate / baseline
##   (so 90 = 10% inhibition; 110 = 10% excitation)
## - Smoothing drives detection (NOT single raw spikes)
## - Persistence gate: excursion must persist >= min_persist_s in window
## - Rapid window biphasic detection (0–10 s) using smoothed + persistence
## - Global extrema computed over post_window using smoothed
##
## Colors (for later plotting):
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

#comp5HT = compile5HT()
df = comp5HT$srtPuff

#unique(df$expCon)

unique(df$assigned_subclass)


# ---- user-side filtering (edit as needed) ----
df_use <- df %>%
  dplyr::filter(
    bath == "none",
    !cell_name %in% c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03", "Q21.26.010.11.13.02"),
    grepl("^5HT", puff),
    expCon == "Standard_Puff"
  ) %>%
  dplyr::mutate(
    assigned_subclass = dplyr::case_when(
      assigned_subclass %in% c("L3c") ~ "L23_IT",
      TRUE ~ assigned_subclass
    )
  ) %>%
  dplyr::filter(assigned_subclass %in% c("L23_IT","L5_ET","L5_IT"))# %>%


df_L6 <- df %>%
  dplyr::filter(
    bath == "none",
    !cell_name %in% c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03"),
    grepl("^5HT", puff),
    expCon == "Standard_Puff"
  ) %>%
  dplyr::mutate(
    assigned_subclass = dplyr::case_when(
      assigned_subclass %in% c("L3c") ~ "L23_IT",
      TRUE ~ assigned_subclass
    )
  ) %>%
  dplyr::filter(grepl("^L6", assigned_subclass))

# =========================================================
# 1) Prep raw df: baseline and pct_of_baseline
# =========================================================
prep_puff_df <- function(df,
                         baseline_window = c(-5, 0)) {

  stopifnot(all(c("cell_name","time","instRate") %in% names(df)))

  df0 <- df %>%
    dplyr::mutate(
      time = as.numeric(time),
      instRate = as.numeric(instRate)
    )

  base_tbl <- df0 %>%
    dplyr::filter(time >= baseline_window[1], time < baseline_window[2]) %>%
    dplyr::group_by(cell_name) %>%
    dplyr::summarise(
      baseline_hz = mean(instRate, na.rm = TRUE),
      n_baseline_pts = sum(is.finite(instRate)),
      .groups = "drop"
    )

  df1 <- df0 %>%
    dplyr::left_join(base_tbl, by = "cell_name") %>%
    dplyr::mutate(
      pct_of_baseline = dplyr::if_else(
        is.finite(baseline_hz) & baseline_hz > 0,
        100 * instRate / baseline_hz,
        NA_real_
      )
    )

  df1
}

# =========================================================
# 2) Smoothing helper
# =========================================================
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

# =========================================================
# 3) Persistence gate for excursions (IMPORTANT)
#    Requires a continuous run beyond threshold for >= min_persist_s
# =========================================================
has_persistent_excursion <- function(d, window, direction = c("exc","inh"),
                                     thr = 5,
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

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  true_runs <- which(r$values)
  if (length(true_runs) == 0) return(FALSE)

  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]
    dur <- d0$time[i1] - d0$time[i0]
    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }

  FALSE
}

time_to_return_baseline <- function(d,
                                    t0,
                                    window = c(0, 20),
                                    col = "pct_smooth",
                                    tol = 3,             # +/- % of baseline (100)
                                    min_persist_s = 0.5) {
  d0 <- d %>%
    dplyr::filter(time >= window[1], time <= window[2]) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(time))

  if (nrow(d0) < 2) return(NA_real_)

  # only look after the peak/trough time (t0)
  d1 <- d0 %>% dplyr::filter(time >= t0)
  if (nrow(d1) < 2) return(NA_real_)

  in_band <- abs(d1[[col]] - 100) <= tol
  if (!any(in_band)) return(NA_real_)

  r <- rle(in_band)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  true_runs <- which(r$values)
  if (!length(true_runs)) return(NA_real_)

  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]
    dur <- d1$time[i1] - d1$time[i0]
    if (is.finite(dur) && dur >= min_persist_s) {
      return(d1$time[i0])   # first time it returned and stayed in band
    }
  }

  NA_real_
}

flag_low_change <- function(per_cell,
                            puff_mag_thr = 5,        # e.g. <5% change at puff
                            max_return_latency = 2,  # return within 2s after peak
                            require_return = TRUE) {

  per_cell %>%
    dplyr::mutate(
      puff_mag_dir = dplyr::case_when(
        direction == "excitation" ~ puff_exc_mag,
        direction == "inhibition" ~ puff_inh_mag,
        TRUE ~ NA_real_
      ),

      low_change_fast_return = dplyr::case_when(
        direction %in% c("excitation","inhibition") ~
          (is.finite(puff_mag_dir) & puff_mag_dir < puff_mag_thr) &
          (if (require_return) is.finite(return_latency_s) else TRUE) &
          (is.finite(return_latency_s) & return_latency_s <= max_return_latency),
        TRUE ~ FALSE
      )
    )
}

early_window <- c(0, 12)
late_window  <- c(12, 56)





# =========================================================
# 4) Feature extraction per cell (SMOOTH + PERSISTENCE)
# =========================================================
extract_bucket_features_one <- function(df_cell,
                                        rapid_window = c(0, 10),
                                        mid_window   = c(10, 40),
                                        late_slope_window = c(36, 56),  # last 20s
                                        late_mean_window  = c(46, 56),  # last 10s
                                        post_window = c(0, 56),
                                        min_excursion_pct = 10,
                                        k_smooth = 11,
                                        smooth_method = c("rollmedian","rollmean","loess"),
                                        loess_span = 0.15,
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

  # --- extrema computed from SMOOTH ---
  rapid_min <- get_extreme(df_cell, rapid_window, "min", col = "pct_smooth")
  rapid_max <- get_extreme(df_cell, rapid_window, "max", col = "pct_smooth")
  abs_min   <- get_extreme(df_cell, post_window,  "min", col = "pct_smooth")
  abs_max   <- get_extreme(df_cell, post_window,  "max", col = "pct_smooth")

  # --- persistence-based excursion flags (SMOOTH) ---
  rapid_has_exc <- has_persistent_excursion(df_cell, rapid_window, "exc",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)
  rapid_has_inh <- has_persistent_excursion(df_cell, rapid_window, "inh",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)

  abs_has_exc <- has_persistent_excursion(df_cell, post_window, "exc",
                                          thr = min_excursion_pct,
                                          col = "pct_smooth",
                                          min_persist_s = min_persist_s)
  abs_has_inh <- has_persistent_excursion(df_cell, post_window, "inh",
                                          thr = min_excursion_pct,
                                          col = "pct_smooth",
                                          min_persist_s = min_persist_s)

  # --- slow windows (for delayed 5HT effects) ---
  early_window <- c(0, 12)
  late_window  <- c(12, 56)

  # --- persistence flags in early/late windows ---
  early_has_exc <- has_persistent_excursion(df_cell, early_window, "exc",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)
  early_has_inh <- has_persistent_excursion(df_cell, early_window, "inh",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)

  late_has_exc  <- has_persistent_excursion(df_cell, late_window, "exc",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)
  late_has_inh  <- has_persistent_excursion(df_cell, late_window, "inh",
                                            thr = min_excursion_pct,
                                            col = "pct_smooth",
                                            min_persist_s = min_persist_s)

  # --- extrema in early/late windows ---
  early_min <- get_extreme(df_cell, early_window, "min", col = "pct_smooth")
  early_max <- get_extreme(df_cell, early_window, "max", col = "pct_smooth")
  late_min  <- get_extreme(df_cell, late_window,  "min", col = "pct_smooth")
  late_max  <- get_extreme(df_cell, late_window,  "max", col = "pct_smooth")

  # --- magnitudes in late window (to avoid tiny rebound being called biphasic) ---
  late_exc_mag <- if (is.finite(late_max$value)) late_max$value - 100 else NA_real_
  late_inh_mag <- if (is.finite(late_min$value)) 100 - late_min$value else NA_real_

  late_mag_thr <- min_excursion_pct  # or set separate e.g. 8–12

  # --- slow pattern label (ONE definition) ---
  slow_pattern <- dplyr::case_when(
    early_has_inh & late_has_exc & is.finite(late_exc_mag) & late_exc_mag >= late_mag_thr ~
      "slow_biphasic_inh_to_exc",

    early_has_exc & late_has_inh & is.finite(late_inh_mag) & late_inh_mag >= late_mag_thr ~
      "slow_biphasic_exc_to_inh",

    early_has_exc & !early_has_inh & late_has_exc ~
      "slow_excitation_sustained_or_delayed",

    early_has_inh & !early_has_exc & late_has_inh ~
      "slow_inhibition_sustained_or_delayed",

    TRUE ~ "slow_no_clear_change"
  )


  # --- gate absolute extrema ---
  max_is_real <- is.finite(abs_max$value) && abs_max$value >= (100 + min_excursion_pct) && abs_has_exc
  min_is_real <- is.finite(abs_min$value) && abs_min$value <= (100 - min_excursion_pct) && abs_has_inh

  abs_max_use <- if (max_is_real) abs_max$value else NA_real_
  abs_min_use <- if (min_is_real) abs_min$value else NA_real_

  exc_mag <- if (is.finite(abs_max_use)) abs_max_use - 100 else NA_real_
  inh_mag <- if (is.finite(abs_min_use)) 100 - abs_min_use else NA_real_

  # late_exc_mag <- if (is.finite(late_max$value)) late_max$value - 100 else NA_real_
  # late_inh_mag <- if (is.finite(late_min$value)) 100 - late_min$value else NA_real_


  late_mag_thr <- min_excursion_pct  # or set separate, e.g. 8 or 10

  direction <- dplyr::case_when(
    is.finite(exc_mag) & is.finite(inh_mag) ~ ifelse(exc_mag >= inh_mag, "excitation", "inhibition"),
    is.finite(exc_mag) & !is.finite(inh_mag) ~ "excitation",
    !is.finite(exc_mag) & is.finite(inh_mag) ~ "inhibition",
    TRUE ~ NA_character_
  )

  direction_slow_aware <- dplyr::case_when(
    slow_pattern == "slow_biphasic_inh_to_exc" ~ "biphasic",
    slow_pattern == "slow_biphasic_exc_to_inh" ~ "biphasic",
    TRUE ~ direction
  )

  response_peak <- dplyr::case_when(
    direction == "excitation" ~ exc_mag,
    direction == "inhibition" ~ inh_mag,
    TRUE ~ NA_real_
  )

  # --- rapid pattern (BIPHASIC uses PERSISTENCE) ---
  rapid_pattern <- dplyr::case_when(
    rapid_has_exc & rapid_has_inh ~
      ifelse(rapid_max$time < rapid_min$time, "rapid_biphasic_exc_to_inh", "rapid_biphasic_inh_to_exc"),
    rapid_has_exc ~ "rapid_excitation_only",
    rapid_has_inh ~ "rapid_inhibition_only",
    TRUE ~ "rapid_no_clear_change"
  )

  # Mid/late summaries from smooth
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

  # ---- low-change logic: small excursion + quick return ----
  # define a tighter "puff window" if you want (e.g., 0–6s)
  puff_window <- c(0, 6)

  puff_min <- get_extreme(df_cell, puff_window, "min", col = "pct_smooth")
  puff_max <- get_extreme(df_cell, puff_window, "max", col = "pct_smooth")

  puff_exc_mag <- if (is.finite(puff_max$value)) puff_max$value - 100 else NA_real_
  puff_inh_mag <- if (is.finite(puff_min$value)) 100 - puff_min$value else NA_real_

  # choose the relevant peak time based on direction (fallback to rapid extrema time)
  # choose peak time by whichever puff excursion is larger (independent of direction)
  peak_time <- dplyr::case_when(
    is.finite(puff_exc_mag) & is.finite(puff_inh_mag) ~
      ifelse(puff_exc_mag >= puff_inh_mag, puff_max$time, puff_min$time),
    is.finite(puff_exc_mag) & !is.finite(puff_inh_mag) ~ puff_max$time,
    !is.finite(puff_exc_mag) & is.finite(puff_inh_mag) ~ puff_min$time,
    TRUE ~ NA_real_
  )

  # fallback to rapid extrema times if puff window isn't informative
  if (!is.finite(peak_time)) {
    peak_time <- dplyr::case_when(
      is.finite(rapid_max$time) & is.finite(rapid_min$time) ~
        ifelse((rapid_max$value - 100) >= (100 - rapid_min$value), rapid_max$time, rapid_min$time),
      is.finite(rapid_max$time) ~ rapid_max$time,
      is.finite(rapid_min$time) ~ rapid_min$time,
      TRUE ~ NA_real_
    )
  }


  return_tol <- 4
  return_window <- c(0, 56)
  return_persist_s <- 3

  t_return <- if (is.finite(peak_time)) {
    time_to_return_baseline(
      df_cell,
      t0 = peak_time,
      window = return_window,
      col = "pct_smooth",
      tol = return_tol,
      min_persist_s = return_persist_s
    )
  } else NA_real_

  # return_latency <- if (is.finite(t_return) && is.finite(peak_time)) t_return - peak_time else NA_real_
  return_latency <- if (is.finite(t_return) && is.finite(peak_time)) t_return - peak_time else NA_real_

  tibble::tibble(
    cell_name = df_cell$cell_name[1],
    assigned_subclass = if ("assigned_subclass" %in% names(df_cell)) df_cell$assigned_subclass[1] else NA_character_,

    baseline_hz = df_cell$baseline_hz[1],
    n_baseline_pts = df_cell$n_baseline_pts[1],

    rapid_min_time = rapid_min$time,
    rapid_min_pct  = rapid_min$value,
    rapid_max_time = rapid_max$time,
    rapid_max_pct  = rapid_max$value,

    abs_min_time = abs_min$time,
    abs_min_pct  = abs_min$value,
    abs_max_time = abs_max$time,
    abs_max_pct  = abs_max$value,

    rapid_has_exc_persist = rapid_has_exc,
    rapid_has_inh_persist = rapid_has_inh,
    abs_has_exc_persist = abs_has_exc,
    abs_has_inh_persist = abs_has_inh,

    min_is_real = min_is_real,
    max_is_real = max_is_real,

    early_has_exc_persist = early_has_exc,
    early_has_inh_persist = early_has_inh,
    late_has_exc_persist  = late_has_exc,
    late_has_inh_persist  = late_has_inh,

    early_min_time = early_min$time,
    early_min_pct  = early_min$value,
    early_max_time = early_max$time,
    early_max_pct  = early_max$value,

    late_min_time = late_min$time,
    late_min_pct  = late_min$value,
    late_max_time = late_max$time,
    late_max_pct  = late_max$value,

    slow_pattern = slow_pattern,


    exc_mag = exc_mag,
    inh_mag = inh_mag,
    direction = direction,
    response_peak = response_peak,

    rapid_pattern = rapid_pattern,

    mid_mean_pct = mid_mean,
    late_mean_pct = late_mean,
    late_slope_pct_per_s = late_slope,

    k_smooth = k_smooth,
    smooth_method = smooth_method,
    min_persist_s = min_persist_s,
    min_excursion_pct = min_excursion_pct,

    puff_min_time = puff_min$time,
    puff_min_pct  = puff_min$value,
    puff_max_time = puff_max$time,
    puff_max_pct  = puff_max$value,

    puff_exc_mag = puff_exc_mag,
    puff_inh_mag = puff_inh_mag,

    peak_time_for_return = peak_time,
    t_return_baseline = t_return,
    return_latency_s = return_latency,

    return_tol_pct = return_tol
  )
}

extract_bucket_features <- function(df_prepped,
                                    rapid_window = c(0, 10),
                                    mid_window   = c(10, 40),
                                    late_slope_window = c(36, 56),
                                    late_mean_window  = c(46, 56),
                                    post_window = c(0, 56),
                                    min_excursion_pct = 10,
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

# =========================================================
# 5) Add rapid-biphasic labels to per_cell (CLEAN + CONSISTENT)
#    Uses rapid_* from the extractor (already persistence-gated)
# =========================================================
add_biphasic_features <- function(per_cell) {

  if (!"slow_pattern" %in% names(per_cell)) per_cell$slow_pattern <- "slow_no_clear_change"

  per_cell %>%
    dplyr::mutate(
      # rapid biphasic (your existing)
      is_biphasic_rapid = rapid_pattern %in% c("rapid_biphasic_exc_to_inh","rapid_biphasic_inh_to_exc"),

      biphasic_order = dplyr::case_when(
        rapid_pattern == "rapid_biphasic_exc_to_inh" ~ "exc_then_inh",
        rapid_pattern == "rapid_biphasic_inh_to_exc" ~ "inh_then_exc",
        slow_pattern  == "slow_biphasic_exc_to_inh"  ~ "exc_then_inh",
        slow_pattern  == "slow_biphasic_inh_to_exc"  ~ "inh_then_exc",
        TRUE ~ NA_character_
      ),

      # biphasic if EITHER rapid or slow says so
      is_biphasic_any = is_biphasic_rapid | slow_pattern %in% c("slow_biphasic_exc_to_inh","slow_biphasic_inh_to_exc"),

      response_type = dplyr::case_when(
        is_biphasic_any ~ "biphasic",
        direction %in% c("excitation","inhibition") ~ direction,
        TRUE ~ "no_change"
      )
    )
}


# =========================================================
# 6) Build analysis tables from df (THIS IS YOUR "BUILD SCRIPT")
# =========================================================
blockers_map <- df_use %>%
  dplyr::distinct(cell_name, blockers)

# ---- prep with baseline [-5,0) ----
df_prepped <- prep_puff_df(df_use, baseline_window = c(-5, 0))

# ---- per-cell features (SMOOTH+PERSISTENCE) ----
per_cell <- extract_bucket_features(
  df_prepped,
  rapid_window = c(0, 10),
  mid_window = c(10, 40),
  late_slope_window = c(36, 56),
  late_mean_window  = c(46, 56),
  post_window = c(0, 56),
  min_excursion_pct = 5,
  k_smooth = 3,
  smooth_method = "rollmean",
  min_persist_s = 0.5
)

# ---- add biphasic labels ----
per_cell2 <- per_cell %>%
  add_biphasic_features() %>%
  flag_low_change(puff_mag_thr = 10, max_return_latency = 20) %>%
  dplyr::mutate(
    response_type = dplyr::case_when(
      response_type %in% c("excitation","inhibition") & low_change_fast_return ~ "no_change",
      TRUE ~ as.character(response_type)
    ),
    response_type = factor(response_type,
                           levels = c("excitation","inhibition","biphasic","no_change"))
  )

# ---- attach blockers + blockers_bin (from df_use) ----
# (per_cell2 may not contain blockers because extractor groups by cell_name only)
per_cell2 <- per_cell2 %>%
  dplyr::left_join(blockers_map, by = "cell_name") %>%
  dplyr::mutate(
    blockers_bin = dplyr::if_else(is.na(blockers) | blockers == "",
                                  "no_blockers", "has_blockers"),
    blockers_bin = factor(blockers_bin, levels = c("no_blockers","has_blockers"))
  )


table(per_cell2$response_type)

## Manual assignments ###

per_cell2 <- per_cell2 %>%
  dplyr::mutate(
    # manual override for this one problematic cell
    response_type = dplyr::if_else(
      cell_name == "Q21.26.017.1A.18.02",
      "inhibition",
      as.character(response_type)
    ),

    is_biphasic_any = dplyr::if_else(
      cell_name == "Q21.26.017.1A.18.02",
      FALSE,
      is_biphasic_any
    ),

    is_biphasic_rapid = dplyr::if_else(
      cell_name == "Q21.26.017.1A.18.02",
      FALSE,
      is_biphasic_rapid
    ),

    biphasic_order = dplyr::if_else(
      cell_name == "Q21.26.017.1A.18.02",
      NA_character_,
      biphasic_order
    )
  ) %>%
  dplyr::mutate(
    response_type = factor(response_type,
                           levels = c("excitation","inhibition","biphasic","no_change")
    )
  )

cells_force_biphasic_inh_then_exc <- c(
  "QN25.26.014.11.12.01",
  "QN26.26.004.20.01.01"
)

per_cell2 <- per_cell2 %>%
  dplyr::mutate(
    response_type = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      "biphasic",
      as.character(response_type)
    ),
    is_biphasic_any = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      TRUE,
      is_biphasic_any
    ),
    is_biphasic_rapid = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      TRUE,
      is_biphasic_rapid
    ),
    biphasic_order = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      "inh_then_exc",
      biphasic_order
    ),
    slow_pattern = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      "slow_biphasic_inh_to_exc",
      slow_pattern
    ),
    rapid_pattern = dplyr::if_else(
      cell_name %in% cells_force_biphasic_inh_then_exc,
      "rapid_biphasic_inh_to_exc",
      rapid_pattern
    )
  ) %>%
  dplyr::mutate(
    response_type = factor(response_type,
                           levels = c("excitation","inhibition","biphasic","no_change"))
  )


# ---- sanity checks you were using ----
per_cell2 %>%
  dplyr::filter(response_type == "excitation") %>%
  dplyr::mutate(
    has_any_inh = is.finite(abs_min_pct) & abs_min_pct <= (100 - 10)
  ) %>%
  dplyr::count(has_any_inh)

per_cell2 %>%
  dplyr::filter(response_type == "excitation",
                is.finite(abs_min_pct),
                abs_min_pct <= 90) %>%
  dplyr::select(cell_name, abs_min_pct, abs_min_time,
                abs_max_pct, abs_max_time, biphasic_order)

# OUTPUTS:
# - df_prepped  : timepoint-level table with baseline_hz + pct_of_baseline
# - per_cell2   : per-cell features + biphasic labels + blockers_bin
