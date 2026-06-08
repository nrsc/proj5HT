# ============================================================
# response_metrics5HT.R — Robust response metrics for 5-HT puffs
# ============================================================
#
# Designed to complement the %-of-baseline classifier in classify5HT.R
# for cells where the baseline is bursty (high CV) or has not fully
# stabilized (drift / adaptation).  Each metric uses a different
# normalisation strategy so that you can cross-check whether a response
# is real, even when the simple percent-change view is dominated by
# baseline noise.
#
# Metrics returned per cell:
#
#   Baseline diagnostics
#     baseline_hz                mean rate in baseline window (Hz)
#     baseline_sd                SD of rate in baseline window (Hz)
#     baseline_cv                CV of rate in baseline window
#     baseline_slope_hz_s        linear drift in baseline (Hz / s)
#
#   Noise-normalised effect size
#     z_peak                     max |z| of response in baseline-SD units
#     z_mean                     mean z of response in baseline-SD units
#     cohens_d                   pooled-SD standardized mean difference
#
#   Drift-corrected magnitude
#     pct_detrended_peak         peak % change vs. extrapolated baseline
#     pct_detrended_mean         mean % change vs. extrapolated baseline
#
#   Changepoint onset
#     cusum_onset_s              time of first CUSUM threshold crossing
#     cusum_direction            "pos" / "neg" / NA
#
#   Regularity change
#     delta_cv_rate              CV(response) - CV(baseline)
#     delta_fano                 Fano(response) - Fano(baseline)
#
#   Distribution shift
#     ks_stat, ks_p              two-sample KS on rate samples
#     wasserstein                1-D earth-mover's distance (Hz)
#
# ============================================================


# ------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------

#' Linearly detrend a baseline and extrapolate through a response window
#'
#' Fits \code{y ~ time} on the baseline window via robust least-squares
#' (ordinary \code{lm}; falls back to mean if too few points), then returns
#' predicted values at the supplied response times.
#'
#' @param t_base,y_base baseline time / rate vectors
#' @param t_resp response times to predict at
#' @return list(predicted = numeric, slope = numeric, intercept = numeric)
#' @keywords internal
.extrapolate_baseline <- function(t_base, y_base, t_resp) {
  ok <- is.finite(t_base) & is.finite(y_base)
  t_base <- t_base[ok]; y_base <- y_base[ok]

  if (length(y_base) < 3) {
    mu <- if (length(y_base) > 0) mean(y_base) else NA_real_
    return(list(
      predicted = rep(mu, length(t_resp)),
      slope     = NA_real_,
      intercept = mu
    ))
  }

  fit <- stats::lm(y_base ~ t_base)
  cf  <- stats::coef(fit)
  pred <- cf[[1]] + cf[[2]] * t_resp

  list(
    predicted = as.numeric(pred),
    slope     = as.numeric(cf[[2]]),
    intercept = as.numeric(cf[[1]])
  )
}

#' One-sided / two-sided CUSUM onset detector
#'
#' Standardises \code{y} by baseline mean/SD then accumulates
#' \eqn{S^{+}_i = \max(0, S^{+}_{i-1} + z_i - k)} and the symmetric
#' negative counterpart.  Returns the first time either |S| crosses
#' \code{h}, with the sign of the crossing.
#'
#' @param time numeric time vector inside the response window
#' @param y    numeric rate samples (same length)
#' @param mu   baseline mean
#' @param sd   baseline SD
#' @param k    slack (in SD units); ignores deviations < k
#' @param h    decision threshold (in SD units)
#' @return list(onset = numeric scalar, direction = "pos"/"neg"/NA)
#' @keywords internal
.cusum_onset <- function(time, y, mu, sd, k = 0.5, h = 5) {
  ok <- is.finite(time) & is.finite(y)
  time <- time[ok]; y <- y[ok]

  if (length(y) == 0 || !is.finite(mu) || !is.finite(sd) || sd <= 0) {
    return(list(onset = NA_real_, direction = NA_character_))
  }

  z <- (y - mu) / sd
  s_pos <- numeric(length(z))
  s_neg <- numeric(length(z))

  for (i in seq_along(z)) {
    prev_p <- if (i == 1) 0 else s_pos[i - 1]
    prev_n <- if (i == 1) 0 else s_neg[i - 1]
    s_pos[i] <- max(0, prev_p + z[i] - k)
    s_neg[i] <- min(0, prev_n + z[i] + k)
  }

  hit_pos <- which(s_pos >=  h)[1]
  hit_neg <- which(s_neg <= -h)[1]

  cand <- c(pos = hit_pos, neg = hit_neg)
  cand <- cand[is.finite(cand)]
  if (length(cand) == 0) {
    return(list(onset = NA_real_, direction = NA_character_))
  }

  first <- which.min(cand)
  list(
    onset     = time[cand[first]],
    direction = names(cand)[first]
  )
}

#' 1-D Wasserstein (earth-mover's) distance between two samples
#'
#' Uses the sorted-quantile form
#' \eqn{W_1(a, b) = \int_0^1 |F_a^{-1}(u) - F_b^{-1}(u)| du}
#' approximated on a common quantile grid.  No new dependencies.
#'
#' @param a,b numeric vectors
#' @param n_grid quantile grid resolution
#' @return numeric scalar (same units as inputs)
#' @keywords internal
.wasserstein1d <- function(a, b, n_grid = 1024) {
  a <- a[is.finite(a)]; b <- b[is.finite(b)]
  if (length(a) == 0 || length(b) == 0) return(NA_real_)

  probs <- seq(0, 1, length.out = n_grid + 2)[-c(1, n_grid + 2)]
  qa <- stats::quantile(a, probs = probs, names = FALSE, na.rm = TRUE)
  qb <- stats::quantile(b, probs = probs, names = FALSE, na.rm = TRUE)
  mean(abs(qa - qb))
}

#' Fano factor of a rate sample (variance / mean), guarded
#' @keywords internal
.fano <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mu <- mean(x)
  if (!is.finite(mu) || mu <= 0) return(NA_real_)
  stats::var(x) / mu
}

#' Coefficient of variation (SD / |mean|), guarded
#' @keywords internal
.cv <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  mu <- mean(x)
  if (!is.finite(mu) || mu == 0) return(NA_real_)
  stats::sd(x) / abs(mu)
}


# ------------------------------------------------------------
# Public per-cell function
# ------------------------------------------------------------

#' Robust response metrics for a single cell's puff trace
#'
#' Computes a panel of complementary response metrics that are robust to
#' bursty or drifting baselines.  Operates on a single cell's rows from
#' the canonical long-format puff dataframe (see \code{prep_puff_df()}).
#'
#' All time values are expected to be aligned so that TTL onset is at
#' \code{time = 0}.
#'
#' @param df_cell data.frame for a single cell with at least the columns
#'   \code{time} (seconds, TTL-aligned) and \code{instRate} (Hz).
#' @param baseline_window numeric length-2.  Pre-puff window (s) used to
#'   estimate baseline mean, SD, and drift.
#' @param response_window numeric length-2.  Post-puff window (s) over
#'   which the response metrics are computed.
#' @param cusum_k slack in baseline-SD units passed to the CUSUM detector.
#' @param cusum_h decision threshold in baseline-SD units for CUSUM.
#' @param wasserstein_grid quantile-grid resolution for the 1-D
#'   Wasserstein computation.
#'
#' @return a one-row \code{tibble} with the metric panel described in
#'   the file header (baseline diagnostics + z / Cohen's d / detrended
#'   percent / CUSUM onset / ΔCV / ΔFano / KS / Wasserstein).
#'
#' @examples
#' \dontrun{
#'   srt <- load5HT("H23.30.001.01")
#'   df  <- srt$dfs$asp$by_protocol[[1]]$spike_puff_output
#'   response_metrics5HT(df)
#' }
#' @export
response_metrics5HT <- function(df_cell,
                                baseline_window = c(-5, 0),
                                response_window = c(0, 10),
                                cusum_k = 0.5,
                                cusum_h = 5,
                                wasserstein_grid = 1024) {

  assert_has_cols(df_cell, c("time", "instRate"))

  df_cell <- df_cell[order(df_cell$time), , drop = FALSE]
  t  <- as.numeric(df_cell$time)
  y  <- as.numeric(df_cell$instRate)

  base_idx <- is.finite(t) & is.finite(y) &
              t >= baseline_window[1] & t <  baseline_window[2]
  resp_idx <- is.finite(t) & is.finite(y) &
              t >= response_window[1] & t <= response_window[2]

  t_base <- t[base_idx]; y_base <- y[base_idx]
  t_resp <- t[resp_idx]; y_resp <- y[resp_idx]

  na_row <- tibble::tibble(
    baseline_hz          = NA_real_,
    baseline_sd          = NA_real_,
    baseline_cv          = NA_real_,
    baseline_slope_hz_s  = NA_real_,
    z_peak               = NA_real_,
    z_mean               = NA_real_,
    cohens_d             = NA_real_,
    pct_detrended_peak   = NA_real_,
    pct_detrended_mean   = NA_real_,
    cusum_onset_s        = NA_real_,
    cusum_direction      = NA_character_,
    delta_cv_rate        = NA_real_,
    delta_fano           = NA_real_,
    ks_stat              = NA_real_,
    ks_p                 = NA_real_,
    wasserstein          = NA_real_
  )

  if (length(y_base) < 3 || length(y_resp) < 2) {
    return(na_row)
  }

  # ---- baseline diagnostics ----
  bl_mean  <- mean(y_base)
  bl_sd    <- stats::sd(y_base)
  bl_cv    <- if (is.finite(bl_mean) && bl_mean > 0) bl_sd / bl_mean else NA_real_
  bl_drift <- tryCatch(
    stats::coef(stats::lm(y_base ~ t_base))[["t_base"]],
    error = function(e) NA_real_
  )

  # ---- noise-normalised effect size ----
  if (is.finite(bl_sd) && bl_sd > 0) {
    z_resp <- (y_resp - bl_mean) / bl_sd
    z_peak <- z_resp[which.max(abs(z_resp))]
    z_mean <- mean(z_resp)
  } else {
    z_peak <- NA_real_; z_mean <- NA_real_
  }

  # Cohen's d with pooled SD
  pooled_sd <- sqrt(
    ((length(y_base) - 1) * stats::var(y_base) +
     (length(y_resp) - 1) * stats::var(y_resp)) /
    (length(y_base) + length(y_resp) - 2)
  )
  cohens_d <- if (is.finite(pooled_sd) && pooled_sd > 0) {
    (mean(y_resp) - bl_mean) / pooled_sd
  } else NA_real_

  # ---- drift-corrected percent change ----
  bl_fit <- .extrapolate_baseline(t_base, y_base, t_resp)
  pred   <- bl_fit$predicted
  pct_dt <- 100 * (y_resp - pred) /
            ifelse(is.finite(pred) & abs(pred) > 1e-9, pred, NA_real_)
  pct_dt_peak <- if (any(is.finite(pct_dt))) pct_dt[which.max(abs(pct_dt))] else NA_real_
  pct_dt_mean <- if (any(is.finite(pct_dt))) mean(pct_dt, na.rm = TRUE) else NA_real_

  # ---- CUSUM onset ----
  cu <- .cusum_onset(t_resp, y_resp, mu = bl_mean, sd = bl_sd,
                     k = cusum_k, h = cusum_h)

  # ---- regularity change ----
  d_cv   <- .cv(y_resp)   - .cv(y_base)
  d_fano <- .fano(y_resp) - .fano(y_base)

  # ---- distribution shift ----
  ks <- tryCatch(
    suppressWarnings(stats::ks.test(y_base, y_resp)),
    error = function(e) NULL
  )
  ks_stat <- if (!is.null(ks)) as.numeric(ks$statistic) else NA_real_
  ks_p    <- if (!is.null(ks)) as.numeric(ks$p.value)   else NA_real_

  wass <- .wasserstein1d(y_base, y_resp, n_grid = wasserstein_grid)

  tibble::tibble(
    baseline_hz          = bl_mean,
    baseline_sd          = bl_sd,
    baseline_cv          = bl_cv,
    baseline_slope_hz_s  = bl_drift,
    z_peak               = z_peak,
    z_mean               = z_mean,
    cohens_d             = cohens_d,
    pct_detrended_peak   = pct_dt_peak,
    pct_detrended_mean   = pct_dt_mean,
    cusum_onset_s        = cu$onset,
    cusum_direction      = cu$direction,
    delta_cv_rate        = d_cv,
    delta_fano           = d_fano,
    ks_stat              = ks_stat,
    ks_p                 = ks_p,
    wasserstein          = wass
  )
}


# ------------------------------------------------------------
# Public group-wise wrapper
# ------------------------------------------------------------

#' Compute response metrics for every cell in a prepared puff dataframe
#'
#' Groups \code{df_prepped} (output of \code{prep_puff_df()}) by
#' \code{cell_name} and applies \code{response_metrics5HT()} to each
#' cell.  The returned per-cell tibble can be left-joined onto the
#' output of \code{build_per_cell_classifier()} to enrich the existing
#' classifier with robust-baseline metrics.
#'
#' @param df_prepped data.frame from \code{prep_puff_df()} (long format,
#'   one row per timepoint per cell).
#' @param baseline_window,response_window,cusum_k,cusum_h,wasserstein_grid
#'   passed through to \code{response_metrics5HT()}.
#'
#' @return tibble with one row per \code{cell_name} and the metric panel.
#' @examples
#' \dontrun{
#'   out      <- classify_puff_response(comp5HT$srtPuff)
#'   metrics  <- build_response_metrics_table(out$df_prepped)
#'   enriched <- dplyr::left_join(out$per_cell, metrics, by = "cell_name")
#' }
#' @export
build_response_metrics_table <- function(df_prepped,
                                         baseline_window = c(-5, 0),
                                         response_window = c(0, 10),
                                         cusum_k = 0.5,
                                         cusum_h = 5,
                                         wasserstein_grid = 1024) {

  assert_has_cols(df_prepped, c("cell_name", "time", "instRate"))

  df_prepped %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(~ response_metrics5HT(
      .x,
      baseline_window  = baseline_window,
      response_window  = response_window,
      cusum_k          = cusum_k,
      cusum_h          = cusum_h,
      wasserstein_grid = wasserstein_grid
    )) %>%
    dplyr::ungroup()
}
