smooth_trace <- function(x, k = 7) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

# Helper: run-length segments of TRUE/FALSE and return indices of TRUE runs
.true_runs <- function(flag) {
  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  data.frame(value = r$values, start = starts, end = ends, len = r$lengths)
}

find_peaks_relative_to_baseline <- function(df_cell,
                                            baseline_window = c(-5, 0),
                                            response_window = c(0, 120),
                                            ignore_after_puff = 0.5,
                                            k_peak = 7,     # light smoothing for peak time
                                            k_sustain = 21, # heavier smoothing for sustained magnitude
                                            sustain_sec = 2,
                                            min_baseline_hz = 0.2) {

  df_cell <- df_cell %>%
    dplyr::arrange(time)

  # Baseline from instRate (Hz)
  base <- df_cell %>%
    dplyr::filter(time >= baseline_window[1], time <= baseline_window[2]) %>%
    dplyr::summarise(b = mean(instRate, na.rm = TRUE)) %>%
    dplyr::pull(b)

  if (is.na(base) || base < min_baseline_hz) {
    return(list(valid = FALSE, reason = "baseline too low/NA", baseline_hz = base))
  }

  d <- df_cell %>%
    dplyr::mutate(
      pc_raw   = 100 * instRate / base,
      pc_peak  = smooth_trace(pc_raw, k = k_peak),
      pc_sust  = smooth_trace(pc_raw, k = k_sustain)
    )

  d_resp <- d %>%
    dplyr::filter(
      time >= (response_window[1] + ignore_after_puff),
      time <= response_window[2]
    )

  if (nrow(d_resp) < 5) {
    return(list(valid = FALSE, reason = "too few points in response window", baseline_hz = base))
  }

  # Peak timing = argmin/argmax of the lightly smoothed curve
  i_inh <- which.min(d_resp$pc_peak)
  i_exc <- which.max(d_resp$pc_peak)

  inh_peak <- d_resp[i_inh, c("time", "pc_peak", "pc_raw", "pc_sust", "instRate")]
  exc_peak <- d_resp[i_exc, c("time", "pc_peak", "pc_raw", "pc_sust", "instRate")]

  # ---- Optional: "sustained segment" logic (still reports the actual min/max within that segment) ----
  # Estimate dt to convert sustain_sec -> points
  dt <- stats::median(diff(d_resp$time), na.rm = TRUE)
  n_need <- max(1, ceiling(sustain_sec / dt))

  # Define sustained inhibition/excitation flags on the heavier-smoothed curve
  # (you can tighten these thresholds later)
  inh_flag <- d_resp$pc_sust <= 95   # at least 5% inhibition (relative to baseline)
  exc_flag <- d_resp$pc_sust >= 105  # at least 5% excitation

  inh_runs <- .true_runs(inh_flag) %>% dplyr::filter(value, len >= n_need)
  exc_runs <- .true_runs(exc_flag) %>% dplyr::filter(value, len >= n_need)

  inh_sust_peak <- NULL
  exc_sust_peak <- NULL

  if (nrow(inh_runs) > 0) {
    # pick the run that contains the deepest minimum in pc_peak
    best <- which.min(sapply(seq_len(nrow(inh_runs)), function(j) {
      rr <- inh_runs[j, ]
      min(d_resp$pc_peak[rr$start:rr$end], na.rm = TRUE)
    }))
    rr <- inh_runs[best, ]
    idx <- rr$start:rr$end
    ii <- idx[which.min(d_resp$pc_peak[idx])]
    inh_sust_peak <- d_resp[ii, c("time", "pc_peak", "pc_raw", "pc_sust", "instRate")]
  }

  if (nrow(exc_runs) > 0) {
    best <- which.max(sapply(seq_len(nrow(exc_runs)), function(j) {
      rr <- exc_runs[j, ]
      max(d_resp$pc_peak[rr$start:rr$end], na.rm = TRUE)
    }))
    rr <- exc_runs[best, ]
    idx <- rr$start:rr$end
    ii <- idx[which.max(d_resp$pc_peak[idx])]
    exc_sust_peak <- d_resp[ii, c("time", "pc_peak", "pc_raw", "pc_sust", "instRate")]
  }

  list(
    valid = TRUE,
    baseline_hz = base,
    df = d,
    inh_peak = inh_peak,
    exc_peak = exc_peak,
    inh_sust_peak = inh_sust_peak,
    exc_sust_peak = exc_sust_peak,
    dt = dt,
    n_need = n_need
  )
}

plot_cell_qc_peak <- function(df, cell_id,
                              baseline_window = c(-5, 0),
                              response_window = c(0, 120),
                              ignore_after_puff = 0.5,
                              k_peak = 7,
                              k_sustain = 21,
                              sustain_sec = 2,
                              use_sustained_peaks = TRUE) {

  df_cell <- df %>% dplyr::filter(cell_name == cell_id)

  out <- find_peaks_relative_to_baseline(
    df_cell,
    baseline_window = baseline_window,
    response_window = response_window,
    ignore_after_puff = ignore_after_puff,
    k_peak = k_peak,
    k_sustain = k_sustain,
    sustain_sec = sustain_sec
  )

  if (!isTRUE(out$valid)) {
    stop(paste("Invalid:", out$reason))
  }

  d <- out$df

  inh_pt <- if (use_sustained_peaks && !is.null(out$inh_sust_peak)) out$inh_sust_peak else out$inh_peak
  exc_pt <- if (use_sustained_peaks && !is.null(out$exc_sust_peak)) out$exc_sust_peak else out$exc_peak

  ggplot(d, aes(time)) +
    annotate("rect",
             xmin = baseline_window[1], xmax = baseline_window[2],
             ymin = -Inf, ymax = Inf, alpha = 0.08, fill = "grey60") +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    geom_line(aes(y = pc_raw), color = "grey75", linewidth = 0.4) +
    geom_line(aes(y = pc_sust), color = "grey35", linewidth = 1.0, alpha = 0.5) +
    geom_line(aes(y = pc_peak), color = "black", linewidth = 0.9) +
    geom_point(data = exc_pt, aes(y = pc_peak), color = "#B2182B", size = 3) +
    geom_point(data = inh_pt, aes(y = pc_peak), color = "#2166AC", size = 3) +
    coord_cartesian(xlim = c(min(d$time, na.rm = TRUE), response_window[2])) +
    labs(
      title = paste0(cell_id, " | baseline=", round(out$baseline_hz, 2), " Hz"),
      subtitle = paste0(
        "inh peak=", round(inh_pt$pc_peak, 1), "% at t=", round(inh_pt$time, 2), "s;  ",
        "exc peak=", round(exc_pt$pc_peak, 1), "% at t=", round(exc_pt$time, 2), "s;  ",
        "dt≈", round(out$dt, 3), "s; sustain=", out$n_need, " pts"
      ),
      x = "Time (s)",
      y = "% of baseline (100 = baseline)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9)
    )
}


plot_cell_qc_peak(df, "QF25.26.024.19.07.04",
                  k_peak = 5, k_sustain = 21,
                  sustain_sec = 2,
                  use_sustained_peaks = TRUE)


library(patchwork)

example_cells <- per_cell_plot %>%
  dplyr::filter(valid) %>%
  dplyr::group_by(assigned_subclass, direction) %>%
  dplyr::slice_sample(n = 1) %>%
  dplyr::ungroup()

plots <- lapply(example_cells$cell_name, function(cid) {
  plot_cell_qc(df, per_cell_plot, cid)
})

wrap_plots(plots, ncol = 3)


qc_cells <- df %>%
  dplyr::distinct(cell_name) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::pull(cell_name)

dir.create("qc_peaks", showWarnings = FALSE)

for (cid in qc_cells) {
  p <- try(
    plot_cell_qc_peak(df, cid, k_peak = 5, k_sustain = 21, sustain_sec = 2),
    silent = TRUE
  )
  if (!inherits(p, "try-error")) {
    ggsave(file.path("qc_peaks", paste0(cid, ".png")), p, width = 7, height = 5, dpi = 150)
  }
}


