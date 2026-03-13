#' Analyze an RMP protocol set (chunk + pulse-train summaries)
#'
#' Computes:
#' - Per-chunk summary stats over aligned time (`tpRMPbl`, `tpRMP_puff`, `trail`)
#' - Optional pulse-train summaries within the `tpRMP_puff` chunk if `pulse_state`
#'   is present (`pulse_on` and `pulse_off` segment features).
#' - Optional plotting artifacts when `return_level = "full"`.
#'
#' @param tst Data frame with at least columns: `time`, `tpRMP`, `chunk`.
#'   If `pulse_state` exists, pulse-train analysis is performed within `chunk == "tpRMP_puff"`.
#' @param on_s Expected pulse ON duration (seconds).
#' @param off_s Expected pulse OFF duration (seconds).
#' @param early_off_s Early window (seconds) within OFF to compute early peak.
#' @param complete_frac Fraction of expected duration required for last segment to be considered complete.
#' @param downsample_for_plot Logical; downsample returned plot dataframe.
#' @param dvdt_thresh_mV_ms Threshold for fast segments in downsampling (if you use adaptive downsample).
#' @param slow_bin_ms Bin size for slow segments in downsampling.
#' @param fast_keep_every Keep every N points in fast segments.
#' @param max_points Cap on returned plot points.
#' @param thin_slow_every Thin factor for slow segments.
#' @param make_plots Logical; if TRUE and full return, includes ggplot trace.
#' @param return_level `"light"` or `"full"`.
#'
#' @return A list with:
#' \describe{
#'   \item{chunk_summary}{tibble of per-chunk summary stats}
#'   \item{bl_mean}{baseline mean (tpRMPbl/tpRMP_bl) or NA}
#'   \item{trail_mean}{trail mean or NA}
#'   \item{pulse}{NULL or list of pulse summaries}
#'   \item{plot_trace_df}{(full) downsampled df for trace plot}
#'   \item{plot_trace}{(full) ggplot trace}
#' }
#'
#' @export
#' @importFrom dplyr filter arrange group_by summarise n pull transmute mutate
#' @importFrom dplyr ungroup left_join select bind_rows slice_tail slice_head
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_point geom_line facet_grid theme_classic theme labs
#' @importFrom stats sd median quantile
analyze_rmp_protocol_set <- function(tst,
                                     on_s  = 0.100,
                                     off_s = 0.900,
                                     early_off_s = 0.100,
                                     complete_frac = 0.80,
                                     downsample_for_plot = TRUE,
                                     dvdt_thresh_mV_ms = 5,
                                     slow_bin_ms = 10,
                                     fast_keep_every = 1,
                                     max_points = 50000,
                                     thin_slow_every = 2,
                                     make_plots = TRUE,
                                     return_level = c("light","full")) {

  return_level <- match.arg(return_level)

  stopifnot(is.data.frame(tst))
  stopifnot(all(c("time","tpRMP","chunk") %in% names(tst)))

  d_all <- tst %>%
    dplyr::filter(is.finite(time), is.finite(tpRMP)) %>%
    dplyr::arrange(time)

  # -----------------------------
  # 1) Chunk summaries
  # -----------------------------
  chunk_summary <- d_all %>%
    dplyr::group_by(chunk) %>%
    dplyr::summarise(
      t_start = min(time, na.rm = TRUE),
      t_end   = max(time, na.rm = TRUE),
      duration_s = t_end - t_start,
      n = dplyr::n(),
      mean_mV = mean(tpRMP, na.rm = TRUE),
      sd_mV   = stats::sd(tpRMP, na.rm = TRUE),
      median_mV = stats::median(tpRMP, na.rm = TRUE),
      p05 = stats::quantile(tpRMP, 0.05, names = FALSE, na.rm = TRUE),
      p95 = stats::quantile(tpRMP, 0.95, names = FALSE, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(t_start)

  # baseline can be named tpRMPbl OR tpRMP_bl depending on upstream
  bl_mean <- chunk_summary %>%
    dplyr::filter(chunk %in% c("tpRMPbl","tpRMP_bl")) %>%
    dplyr::pull(mean_mV)
  trail_mean <- chunk_summary %>%
    dplyr::filter(chunk == "trail") %>%
    dplyr::pull(mean_mV)

  bl_mean <- if (length(bl_mean) == 0) NA_real_ else bl_mean[1]
  trail_mean <- if (length(trail_mean) == 0) NA_real_ else trail_mean[1]

  # -----------------------------
  # 2) Pulse-train analysis within tpRMP_puff
  # -----------------------------
  pulse_res <- NULL
  if ("pulse_state" %in% names(d_all)) {

    d_puff <- d_all %>%
      dplyr::filter(chunk == "tpRMP_puff") %>%
      dplyr::filter(!is.na(pulse_state), pulse_state %in% c("pulse_on","pulse_off")) %>%
      dplyr::transmute(epoch_t = time, tpRMP = tpRMP, pulse_state = pulse_state)

    if (nrow(d_puff) > 0) {

      d <- d_puff %>%
        dplyr::arrange(epoch_t) %>%
        dplyr::mutate(
          state_change = pulse_state != dplyr::lag(pulse_state, default = first(pulse_state)),
          pulse_run    = cumsum(state_change),
          pulse_on_idx  = dplyr::if_else(pulse_state == "pulse_on",
                                         cumsum(state_change & pulse_state == "pulse_on"),
                                         NA_integer_),
          pulse_off_idx = dplyr::if_else(pulse_state == "pulse_off",
                                         cumsum(state_change & pulse_state == "pulse_off"),
                                         NA_integer_),
          pulse_id = dplyr::case_when(
            pulse_state == "pulse_on"  ~ sprintf("on_%03d",  pulse_on_idx),
            pulse_state == "pulse_off" ~ sprintf("off_%03d", pulse_off_idx)
          )
        ) %>%
        dplyr::group_by(pulse_id) %>%
        dplyr::mutate(t0 = min(epoch_t, na.rm = TRUE),
                      t_rel = epoch_t - t0) %>%
        dplyr::ungroup()

      seg_info <- d %>%
        dplyr::group_by(pulse_id, pulse_state) %>%
        dplyr::summarise(
          t_start = min(epoch_t, na.rm = TRUE),
          t_end   = max(epoch_t, na.rm = TRUE),
          duration = t_end - t_start,
          .groups = "drop"
        ) %>%
        dplyr::arrange(t_start)

      # --- drop trailing incomplete ---
      if (nrow(seg_info) >= 1) {
        last_seg <- seg_info %>% dplyr::slice_tail(n = 1)
        exp_last <- if (last_seg$pulse_state == "pulse_on") on_s else off_s
        last_incomplete <- is.finite(last_seg$duration) &&
          (last_seg$duration < complete_frac * exp_last)

        drop_ids <- character(0)

        if (last_incomplete) {
          drop_ids <- last_seg$pulse_id
          if (last_seg$pulse_state == "pulse_on" && nrow(seg_info) >= 2) {
            drop_ids <- c(
              drop_ids,
              seg_info %>%
                dplyr::slice_tail(n = 2) %>%
                dplyr::slice_head(n = 1) %>%
                dplyr::pull(pulse_id)
            )
          }
        } else {
          if (last_seg$pulse_state == "pulse_on" && nrow(seg_info) >= 2) {
            drop_ids <- seg_info %>% dplyr::slice_tail(n = 2) %>% dplyr::pull(pulse_id)
          }
        }

        if (length(drop_ids) > 0) {
          d <- d %>% dplyr::filter(!pulse_id %in% drop_ids)
          seg_info <- seg_info %>% dplyr::filter(!pulse_id %in% drop_ids)
        }
      }

      pulse_on_summary <- d %>%
        dplyr::filter(pulse_state == "pulse_on") %>%
        dplyr::group_by(pulse_id) %>%
        dplyr::summarise(
          mean_on = mean(tpRMP, na.rm = TRUE),
          peak_neg = min(tpRMP, na.rm = TRUE),
          t_peak   = t_rel[which.min(tpRMP)],
          duration = max(t_rel, na.rm = TRUE),
          delta_mean_minus_peak = mean_on - peak_neg,
          .groups = "drop"
        )

      pulse_tau <- d %>%
        dplyr::filter(pulse_state == "pulse_on") %>%
        dplyr::group_by(pulse_id) %>%
        dplyr::summarise(
          peak_tpRMP = min(tpRMP, na.rm = TRUE),
          steady_tpRMP = mean(tpRMP[t_rel > 0.8 * max(t_rel, na.rm = TRUE)], na.rm = TRUE),
          target = peak_tpRMP + 0.63 * (steady_tpRMP - peak_tpRMP),
          tau = t_rel[which.min(abs(tpRMP - target))],
          .groups = "drop"
        )

      pulse_off_summary <- d %>%
        dplyr::filter(pulse_state == "pulse_off") %>%
        dplyr::group_by(pulse_id) %>%
        dplyr::summarise(
          mean_off = mean(tpRMP, na.rm = TRUE),
          sd_off   = stats::sd(tpRMP, na.rm = TRUE),
          duration = max(t_rel, na.rm = TRUE),
          early_window_s = early_off_s,
          early_peak_off = {
            idx <- which(t_rel <= early_window_s & is.finite(tpRMP))
            if (length(idx) == 0) NA_real_ else max(tpRMP[idx], na.rm = TRUE)
          },
          early_t_peak_off = {
            idx <- which(t_rel <= early_window_s & is.finite(tpRMP))
            if (length(idx) == 0) NA_real_ else t_rel[idx][which.max(tpRMP[idx])]
          },
          .groups = "drop"
        )

      if (return_level == "full") {
        pulse_res <- list(
          seg_info = seg_info,
          pulse_on = pulse_on_summary,
          pulse_tau = pulse_tau,
          pulse_off = pulse_off_summary,
          # plot built by separate exported plotting helper
          feat_long = NULL,
          plot_features = NULL,
          df_puff = d
        )
      } else {
        pulse_res <- list(
          seg_info = seg_info,
          pulse_on = pulse_on_summary,
          pulse_tau = pulse_tau,
          pulse_off = pulse_off_summary
        )
      }
    }
  }

  # -----------------------------
  # 3) plotting artifacts (trace) ONLY if full
  # -----------------------------
  df_plot <- NULL
  p_trace <- NULL

  if (return_level == "full") {
    df_plot <- d_all

    # NOTE: if you want adaptive downsampling, keep that helper in helpers,
    # and call it here. For now we keep this minimal and deterministic.
    if (isTRUE(downsample_for_plot) && nrow(df_plot) > max_points) {
      keep <- unique(round(seq(1, nrow(df_plot), length.out = max_points)))
      df_plot <- df_plot[keep, , drop = FALSE]
    }

    if (isTRUE(make_plots) && !is.null(df_plot)) {
      p_trace <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = tpRMP, color = chunk)) +
        ggplot2::geom_line(linewidth = 0.35, alpha = 0.9) +
        ggplot2::theme_classic() +
        ggplot2::labs(x = "Aligned time (s)", y = "tpRMP (mV)", color = "Chunk") +
        ggplot2::theme(legend.position = "top")
    }
  }

  out <- list(
    chunk_summary = chunk_summary,
    bl_mean = bl_mean,
    trail_mean = trail_mean,
    pulse = pulse_res
  )

  if (return_level == "full") {
    out$plot_trace_df <- df_plot
    out$plot_trace <- p_trace
  }

  out
}

#' Plot protocol summary trace across baseline, puff, and trail
#'
#' Convenience plotting wrapper for the aligned `tst` table produced upstream
#' (your concatenated/shifted protocol set). This plots *all chunks* (baseline,
#' puff, trail). If `pulse_state` exists, puff points are lightly overlaid to
#' show pulse_on vs pulse_off.
#'
#' @param tst Data frame with columns `time`, `tpRMP`, `chunk`, and optionally `pulse_state`.
#' @param overlay_pulse_state Logical; if TRUE and `pulse_state` exists, overlays points on puff chunk.
#' @param alpha_line Alpha for trace line.
#' @param linewidth Line width.
#' @return A ggplot object.
#'
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_line geom_point theme_classic theme labs
plot_rmp_protocol_summary <- function(tst,
                                      overlay_pulse_state = TRUE,
                                      alpha_line = 0.90,
                                      linewidth = 0.35) {
  stopifnot(is.data.frame(tst))
  stopifnot(all(c("time","tpRMP","chunk") %in% names(tst)))

  p <- ggplot2::ggplot(tst, ggplot2::aes(x = time, y = tpRMP, color = chunk)) +
    ggplot2::geom_line(linewidth = linewidth, alpha = alpha_line) +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "Aligned time (s)", y = "tpRMP (mV)", color = "Chunk") +
    ggplot2::theme(legend.position = "top")

  if (isTRUE(overlay_pulse_state) && ("pulse_state" %in% names(tst))) {
    puff <- tst %>% dplyr::filter(chunk == "tpRMP_puff" & !is.na(pulse_state))
    if (nrow(puff) > 0) {
      # overlay as points; pulse_state mapped to shape (keeps chunk colors intact)
      p <- p + ggplot2::geom_point(
        data = puff,
        mapping = ggplot2::aes(shape = pulse_state),
        alpha = 0.35,
        size = 0.6,
        inherit.aes = TRUE
      )
    }
  }

  p
}

#' Plot pulse-train feature summaries for an RMP protocol
#'
#' @param pulse A list produced by `analyze_rmp_protocol_set()$pulse`.
#'   In `return_level="full"` this includes `feat_long` and `plot_features`.
#'   In `return_level="light"` this includes only compact tables.
#' @param prefer_cached Logical; if TRUE and `pulse$plot_features` exists, return it.
#'
#' @return A ggplot object (or NULL if no pulse data).
#' @export
plot_rmp_pulse_features <- function(pulse, prefer_cached = TRUE) {
  if (is.null(pulse)) return(NULL)
  stopifnot(is.list(pulse))

  # If full mode already built the ggplot, just use it
  if (isTRUE(prefer_cached) && !is.null(pulse$plot_features)) {
    return(pulse$plot_features)
  }

  # Otherwise rebuild from what we have.
  # LIGHT mode does NOT include feat_long, so we rebuild a compact long table.
  if (!is.null(pulse$feat_long)) {
    feat_long <- pulse$feat_long
  } else {
    # Expect these tables in both light+full
    if (is.null(pulse$pulse_on) && is.null(pulse$pulse_off)) return(NULL)

    on_long <- NULL
    if (!is.null(pulse$pulse_on) && nrow(pulse$pulse_on) > 0) {
      on_long <- pulse$pulse_on |>
        dplyr::mutate(
          pulse_state = "pulse_on",
          pulse_idx = as.integer(sub("^on_", "", pulse_id))
        ) |>
        dplyr::select(pulse_state, pulse_idx,
                      peak_neg, mean_on, delta_mean_minus_peak, t_peak, duration) |>
        tidyr::pivot_longer(
          cols = c(peak_neg, mean_on, delta_mean_minus_peak, t_peak, duration),
          names_to = "feature", values_to = "value"
        )
    }

    off_long <- NULL
    if (!is.null(pulse$pulse_off) && nrow(pulse$pulse_off) > 0) {
      off_long <- pulse$pulse_off |>
        dplyr::mutate(
          pulse_state = "pulse_off",
          pulse_idx = as.integer(sub("^off_", "", pulse_id))
        ) |>
        dplyr::select(pulse_state, pulse_idx,
                      mean_off, early_peak_off, early_t_peak_off, duration) |>
        tidyr::pivot_longer(
          cols = c(mean_off, early_peak_off, early_t_peak_off, duration),
          names_to = "feature", values_to = "value"
        )
    }

    feat_long <- dplyr::bind_rows(on_long, off_long)
    if (is.null(feat_long) || nrow(feat_long) == 0) return(NULL)
  }

  ggplot2::ggplot(feat_long, ggplot2::aes(x = pulse_idx, y = value, color = pulse_state)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.85) +
    ggplot2::geom_line(ggplot2::aes(group = interaction(pulse_state, feature)),
                       linewidth = 0.4, alpha = 0.6) +
    ggplot2::facet_grid(feature ~ ., scales = "free_y", switch = "y") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::labs(x = "Pulse index", y = NULL, color = "State")
}

#' Plot an RMP protocol object (trace or features)
#'
#' @param rmp_protocol One element of `srt$dfs$rmp$by_protocol`.
#' @param what "trace" or "features".
#'
#' @return A ggplot object (or NULL).
#' @export
plot_rmp_protocol <- function(rmp_protocol, what = c("trace","features")) {
  what <- match.arg(what)
  stopifnot(is.list(rmp_protocol))

  if (what == "features") {
    return(plot_rmp_pulse_features(rmp_protocol$rmp_protocol_analysis$pulse))
  }

  # trace requires trace-level df
  tst <- rmp_protocol$rmp_protocol_analysis$plot_trace_df %||%
    rmp_protocol$rmp_protocol_analysis$tst %||% NULL

  if (is.null(tst)) {
    stop("No trace df stored. Rebuild with store_level='full' or store a thinned tst.")
  }

  plot_rmp_protocol_summary(tst)
}

