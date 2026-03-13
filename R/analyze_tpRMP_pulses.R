analyze_tpRMP_pulses <- function(df,
                                 sr = NULL,
                                 on_s  = 0.100,
                                 off_s = 0.900,
                                 early_off_s = 0.100,
                                 # completeness tolerance (fraction of expected duration)
                                 complete_frac = 0.80,
                                 # downsample settings
                                 downsample_for_plot = TRUE,
                                 dvdt_thresh_mV_ms = 5,
                                 slow_bin_ms = 10,
                                 fast_keep_every = 1,
                                 max_points = 50000,
                                 thin_slow_every = 2) {

  stopifnot(is.data.frame(df))
  stopifnot(all(c("epoch_t", "tpRMP", "pulse_state") %in% names(df)))

  # ---- 1) order + build pulse_id + t_rel ----
  d <- df %>%
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
    dplyr::mutate(
      t0 = min(epoch_t, na.rm = TRUE),
      t_rel = epoch_t - t0
    ) %>%
    dplyr::ungroup()

  # ---- 2) drop trailing incomplete segments based on expected on/off duration ----
  seg_info <- d %>%
    dplyr::group_by(pulse_id, pulse_state) %>%
    dplyr::summarise(
      t_start = min(epoch_t, na.rm = TRUE),
      t_end   = max(epoch_t, na.rm = TRUE),
      duration = t_end - t_start,
      .groups = "drop"
    ) %>%
    dplyr::arrange(t_start)

  if (nrow(seg_info) >= 1) {
    last_seg <- seg_info %>% dplyr::slice_tail(n = 1)

    exp_last <- if (last_seg$pulse_state == "pulse_on") on_s else off_s
    last_incomplete <- is.finite(last_seg$duration) &&
      (last_seg$duration < complete_frac * exp_last)

    drop_ids <- character(0)

    if (last_incomplete) {
      drop_ids <- last_seg$pulse_id

      # If we end on an incomplete ON, also drop the preceding OFF so we "end on complete OFF"
      if (last_seg$pulse_state == "pulse_on" && nrow(seg_info) >= 2) {
        drop_ids <- c(drop_ids, seg_info %>% dplyr::slice_tail(n = 2) %>% dplyr::slice_head(n = 1) %>% dplyr::pull(pulse_id))
      }
    } else {
      # If last segment is complete but it's an ON, you said you want to end on OFF
      if (last_seg$pulse_state == "pulse_on" && nrow(seg_info) >= 2) {
        drop_ids <- seg_info %>% dplyr::slice_tail(n = 2) %>% dplyr::pull(pulse_id)
      }
    }

    if (length(drop_ids) > 0) {
      d <- d %>% dplyr::filter(!pulse_id %in% drop_ids)
      seg_info <- seg_info %>% dplyr::filter(!pulse_id %in% drop_ids)
    }
  }

  # ---- 3) pulse ON features ----
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

  # ---- 4) pulse OFF features (incl early peak in first 100 ms) ----
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

  # ---- 5) tidy long + one ggplot with facets ----
  on_df <- pulse_on_summary %>%
    dplyr::left_join(pulse_tau, by = "pulse_id") %>%
    dplyr::mutate(
      pulse_state = "pulse_on",
      pulse_idx = as.integer(sub("^on_", "", pulse_id))
    ) %>%
    dplyr::select(pulse_state, pulse_idx, peak_neg, mean_on, delta_mean_minus_peak, t_peak, tau)

  off_df <- pulse_off_summary %>%
    dplyr::mutate(
      pulse_state = "pulse_off",
      pulse_idx = as.integer(sub("^off_", "", pulse_id))
    ) %>%
    dplyr::select(pulse_state, pulse_idx, mean_off, early_peak_off, early_t_peak_off)

  feat_long <- dplyr::bind_rows(
    on_df %>%
      tidyr::pivot_longer(
        cols = c(peak_neg, mean_on, delta_mean_minus_peak, t_peak, tau),
        names_to = "feature", values_to = "value"
      ),
    off_df %>%
      tidyr::pivot_longer(
        cols = c(mean_off, early_peak_off, early_t_peak_off),
        names_to = "feature", values_to = "value"
      )
  ) %>%
    dplyr::mutate(
      feature = factor(
        feature,
        levels = c(
          "peak_neg", "mean_on", "delta_mean_minus_peak", "t_peak", "tau",
          "mean_off", "early_peak_off", "early_t_peak_off"
        ),
        labels = c(
          "Pulse ON: peak negative (mV)",
          "Pulse ON: mean (mV)",
          "Pulse ON: mean − peak (mV)",
          "Pulse ON: time to peak (s)",
          "Pulse ON: tau to steady (s)",
          "Pulse OFF: mean (mV)",
          "Pulse OFF: early peak ≤100ms (mV)",
          "Pulse OFF: time to early peak (s)"
        )
      ),
      pulse_state = factor(pulse_state, levels = c("pulse_off", "pulse_on"))
    )

  p_all <- ggplot2::ggplot(feat_long, ggplot2::aes(x = pulse_idx, y = value, color = pulse_state)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.85) +
    ggplot2::geom_line(
      ggplot2::aes(group = interaction(pulse_state, feature)),
      linewidth = 0.4, alpha = 0.6
    ) +
    ggplot2::facet_grid(feature ~ ., scales = "free_y", switch = "y") +
    ggplot2::labs(x = "Pulse number", y = NULL, color = "State") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0),
      legend.position = "top"
    )

  # ---- 6) downsampled trace for lightweight plotting only ----
  df_plot <- NULL
  if (isTRUE(downsample_for_plot)) {
    # pick sr if not provided (try infer from epoch_t spacing)
    if (is.null(sr)) {
      dt_med <- stats::median(diff(d$epoch_t), na.rm = TRUE)
      if (is.finite(dt_med) && dt_med > 0) sr <- 1 / dt_med
    }

    if (!is.null(sr) && is.finite(sr) && sr > 0) {
      ds <- downsample_trace_adaptive(
        v = d$tpRMP,
        sr = sr,
        dvdt_thresh_mV_ms = dvdt_thresh_mV_ms,
        slow_bin_ms = slow_bin_ms,
        fast_keep_every = fast_keep_every,
        max_points = max_points,
        thin_slow_every = thin_slow_every
      )
      # keep state + time in same downsampled rows
      df_plot <- dplyr::tibble(
        epoch_t = d$epoch_t[ds$idx],
        tpRMP = d$tpRMP[ds$idx],
        pulse_state = d$pulse_state[ds$idx]
      )
    } else {
      # fallback: thin every Nth point if sr can't be inferred
      keep <- seq(1, nrow(d), by = 10)
      df_plot <- dplyr::tibble(
        epoch_t = d$epoch_t[keep],
        tpRMP = d$tpRMP[keep],
        pulse_state = d$pulse_state[keep]
      )
    }
  }

  # ---- outputs ----
  list(
    # big df NOT returned (lightweight): you can toggle if needed later
    seg_info = seg_info,
    pulse_on = pulse_on_summary,
    pulse_tau = pulse_tau,
    pulse_off = pulse_off_summary,
    feat_long = feat_long,
    plot_features = p_all,
    plot_trace_df = df_plot
  )
}
