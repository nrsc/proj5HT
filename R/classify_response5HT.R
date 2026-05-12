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

      any_exc <- rapid_exc || ext_exc
      any_inh <- rapid_inh || ext_inh

      # extrema in the combined window — used to order biphasic
      ex_max <- get_extreme_in_window(t, ys, combined_window, "max")
      ex_min <- get_extreme_in_window(t, ys, combined_window, "min")

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
