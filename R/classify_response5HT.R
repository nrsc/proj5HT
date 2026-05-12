#' Classify 5-HT puff responses into 4 categories
#'
#' A focused, self-contained classifier that uses both a "rapid" early
#' window and an "extended" later window to call each cell as one of:
#'
#'   * `"excitation"`        — sustained excursion above baseline only
#'   * `"inhibition"`        — sustained excursion below baseline only
#'   * `"biphasic_exc_inh"`  — both directions, excitation comes first
#'   * `"biphasic_inh_exc"`  — both directions, inhibition comes first
#'   * `"no_change"`         — no sustained excursion meeting threshold
#'
#' Operates directly on a long-format dataframe (one row per cell × time)
#' with at minimum `cell_name`, `time`, `percent_change`. This is the same
#' format returned by `figure_df_build()`, so the classifier can be applied
#' to any of those dfs without a separate prep step.
#'
#' A "sustained excursion" is a contiguous run of smoothed
#' `percent_change` that crosses `100 ± thr_pct` for at least
#' `min_persist_s` seconds inside the given window.
#'
#' In addition to the smoothed test, an unsmoothed **zero-spiking** check
#' is run: any contiguous run of raw `percent_change <= zero_thr_pct` (i.e.
#' the cell stopped firing) lasting at least `min_zero_s` seconds inside
#' a window forces the corresponding `*_inh` flag to TRUE. This guarantees
#' that a brief but unambiguous silent period is never missed by the
#' smoother, so an otherwise excitatory cell with a transient silence at
#' puff onset is correctly classified as biphasic.
#'
#' The biphasic order is decided by the order of the two extrema (max for
#' excitation, min for inhibition) inside the full `combined_window`. The
#' rapid window matters because it gives a stricter detection bar for
#' fast responses, but a slow inhibitory rebound in the extended window
#' still counts toward biphasic.
#'
#' @param df      Long-format dataframe. Must contain `cell_name`, `time`,
#'                `percent_change`. Optional grouping columns (e.g.
#'                `assigned_subclass`, `Species`) are carried through.
#' @param rapid_window     numeric length-2, default `c(0, 10)`.
#' @param extended_window  numeric length-2, default `c(10, 50)`.
#' @param thr_pct          excursion threshold (percent units), default `15`.
#' @param min_persist_s    minimum contiguous duration (s), default `1`.
#' @param smooth_k         odd integer rolling-mean window, default `5`.
#' @param zero_thr_pct     `percent_change` below this is treated as silent
#'                         (default `5`).
#' @param min_zero_s       minimum contiguous silent duration (s) to flag
#'                         inhibition (default `0.5`).
#' @param carry_cols       columns to keep on the per-cell output. Auto-
#'                         detected from `df` if `NULL`.
#'
#' @return tibble with one row per cell containing:
#'   `cell_name`, optional carried columns, `response_type` (factor),
#'   `rapid_exc`, `rapid_inh`, `ext_exc`, `ext_inh` (logical persistence
#'   flags), `exc_peak_time`, `inh_peak_time`, `exc_peak_pct`,
#'   `inh_peak_pct`.
#'
#' @seealso [classify5HT_helpers.R] for the `smooth_vec`,
#'   `has_persistent_excursion`, and `get_extreme_in_window` helpers reused
#'   here.
#' @export
classify_response5HT <- function(df,
                                 rapid_window    = c(0, 10),
                                 extended_window = c(10, 50),
                                 thr_pct         = 15,
                                 min_persist_s   = 1,
                                 smooth_k        = 5,
                                 zero_thr_pct    = 5,
                                 min_zero_s      = 0.5,
                                 carry_cols      = NULL) {

  assert_has_cols(df, c("cell_name", "time", "percent_change"))

  if (is.null(carry_cols)) {
    carry_cols <- intersect(
      c("assigned_subclass", "Species", "bath", "puff", "expCon",
        "blockers", "cluster_Corr", "assigned_depth", "depth_bin"),
      names(df)
    )
  }

  combined_window <- c(min(rapid_window[1], extended_window[1]),
                       max(rapid_window[2], extended_window[2]))

  per_cell <- df %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(function(.x, .y) {
      .x <- .x %>% dplyr::arrange(.data$time)
      t  <- as.numeric(.x$time)
      y  <- as.numeric(.x$percent_change)
      ys <- smooth_vec(y, method = "rollmean", k = smooth_k)

      # persistence flags in each window
      rapid_exc <- has_persistent_excursion(t, ys, rapid_window,    "exc",
                                            baseline = 100,
                                            thr      = thr_pct,
                                            min_persist_s = min_persist_s)
      rapid_inh <- has_persistent_excursion(t, ys, rapid_window,    "inh",
                                            baseline = 100,
                                            thr      = thr_pct,
                                            min_persist_s = min_persist_s)
      ext_exc   <- has_persistent_excursion(t, ys, extended_window, "exc",
                                            baseline = 100,
                                            thr      = thr_pct,
                                            min_persist_s = min_persist_s)
      ext_inh   <- has_persistent_excursion(t, ys, extended_window, "inh",
                                            baseline = 100,
                                            thr      = thr_pct,
                                            min_persist_s = min_persist_s)

      # ----- raw zero-spiking detector (bypasses smoothing) -----
      # A sustained run of raw percent_change <= zero_thr_pct means the cell
      # actually stopped firing. Smoothing can hide such brief silences, so
      # we force the inhibition flag for whichever window the silence lies in.
      rapid_zero <- .has_zero_run(t, y, rapid_window,
                                  zero_thr_pct, min_zero_s)
      ext_zero   <- .has_zero_run(t, y, extended_window,
                                  zero_thr_pct, min_zero_s)
      rapid_inh  <- rapid_inh || rapid_zero
      ext_inh    <- ext_inh   || ext_zero

      any_exc <- rapid_exc || ext_exc
      any_inh <- rapid_inh || ext_inh

      # extrema in the combined window — used to order biphasic
      ex_max <- get_extreme_in_window(t, ys, combined_window, "max")
      ex_min <- get_extreme_in_window(t, ys, combined_window, "min")

      # If the inhibition flag was set by a raw hard-zero rather than by a
      # smoothed dip, the smoothed `ex_min` will not correspond to the true
      # silent moment (the rollmean averages a single zero into its bright
      # neighbours). Use the time of the earliest raw zero/near-zero sample
      # instead, so biphasic ordering reflects the actual silence.
      if (rapid_zero || ext_zero) {
        zwin <- if (rapid_zero) rapid_window else extended_window
        raw_zero_time <- .first_zero_time(t, y, zwin, zero_thr_pct,
                                          hard_zero_pct = 1)
        if (is.finite(raw_zero_time)) ex_min$time <- raw_zero_time
      }

      response_type <- if (any_exc && any_inh) {
        # decide order by which extremum occurred first
        if (is.finite(ex_max$time) && is.finite(ex_min$time)) {
          if (ex_max$time <= ex_min$time) "biphasic_exc_inh" else "biphasic_inh_exc"
        } else if (rapid_exc && !rapid_inh) {
          "biphasic_exc_inh"
        } else if (rapid_inh && !rapid_exc) {
          "biphasic_inh_exc"
        } else {
          # fall back: rapid window flagged both — call by which threshold
          # was crossed first in the rapid window
          .order_by_first_cross(t, ys, rapid_window, thr_pct)
        }
      } else if (any_exc) {
        "excitation"
      } else if (any_inh) {
        "inhibition"
      } else {
        "no_change"
      }

      tibble::tibble(
        rapid_exc       = rapid_exc,
        rapid_inh       = rapid_inh,
        ext_exc         = ext_exc,
        ext_inh         = ext_inh,
        rapid_zero      = rapid_zero,
        ext_zero        = ext_zero,
        exc_peak_time   = ex_max$time,
        exc_peak_pct    = ex_max$value,
        inh_peak_time   = ex_min$time,
        inh_peak_pct    = ex_min$value,
        response_type   = response_type
      )
    }) %>%
    dplyr::ungroup()

  # carry through per-cell metadata
  if (length(carry_cols) > 0) {
    meta <- df %>%
      dplyr::group_by(.data$cell_name) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(carry_cols),
                                     ~ first_non_missing(.x)),
                       .groups = "drop")
    per_cell <- dplyr::left_join(per_cell, meta, by = "cell_name")
  }

  per_cell %>%
    dplyr::mutate(
      response_type = factor(
        .data$response_type,
        levels = c("excitation", "inhibition",
                   "biphasic_exc_inh", "biphasic_inh_exc",
                   "no_change")
      )
    )
}


#' Time of the earliest raw zero / hard-zero sample inside a window
#'
#' Used to anchor biphasic ordering when the inhibition flag was set by
#' [`.has_zero_run`] rather than by a smoothed dip.
#'
#' @keywords internal
.first_zero_time <- function(t, y, window, zero_thr_pct, hard_zero_pct = 1) {
  keep <- is.finite(t) & is.finite(y) &
    t >= window[1] & t <= window[2]
  t0 <- t[keep]; y0 <- y[keep]
  if (length(t0) == 0) return(NA_real_)
  hits <- which(y0 <= hard_zero_pct)
  if (length(hits) == 0) hits <- which(y0 <= zero_thr_pct)
  if (length(hits) == 0) return(NA_real_)
  t0[hits[1]]
}

#' Detect a sustained near-zero run of raw `percent_change` inside a window
#'
#' Bypasses smoothing — used to catch brief silent periods that the rolling
#' mean might wash out.
#'
#' Two paths flip the result to TRUE:
#'   * a contiguous run of `y <= zero_thr_pct` lasting at least `min_zero_s`
#'     (the standard persistence rule), OR
#'   * any single raw sample at or below `hard_zero_pct` (default `1`). A
#'     literal zero / near-zero in the trace is treated as unambiguous
#'     evidence of silence — even one sample's worth — because it can't
#'     plausibly be produced by noise. This catches cells that drop to a
#'     single hard zero at puff onset and then rebound to a high firing
#'     rate before the next sample (a pattern the run-length test misses
#'     on coarse bin grids).
#'
#' @keywords internal
.has_zero_run <- function(t, y, window, zero_thr_pct, min_zero_s,
                          hard_zero_pct = 1) {
  keep <- is.finite(t) & is.finite(y) &
    t >= window[1] & t <= window[2]
  t0 <- t[keep]; y0 <- y[keep]
  if (length(t0) == 0) return(FALSE)

  # 1. unambiguous hard zero — any single sample is enough
  if (any(y0 <= hard_zero_pct, na.rm = TRUE)) return(TRUE)

  if (length(t0) < 2) return(FALSE)

  # 2. standard persistence rule for noisier near-zero dips
  runs <- find_runs(t0, y0 <= zero_thr_pct)
  if (nrow(runs) == 0) return(FALSE)
  any(runs$dur >= min_zero_s, na.rm = TRUE)
}

#' Tie-breaker for biphasic ordering when both extrema are in the rapid window
#' @keywords internal
.order_by_first_cross <- function(t, y, window, thr) {
  keep <- is.finite(t) & is.finite(y) & t >= window[1] & t <= window[2]
  t0 <- t[keep]; y0 <- y[keep]
  if (length(t0) == 0) return("biphasic_exc_inh")
  exc_idx <- which(y0 >= 100 + thr)
  inh_idx <- which(y0 <= 100 - thr)
  t_exc <- if (length(exc_idx)) t0[exc_idx[1]] else Inf
  t_inh <- if (length(inh_idx)) t0[inh_idx[1]] else Inf
  if (t_exc <= t_inh) "biphasic_exc_inh" else "biphasic_inh_exc"
}
