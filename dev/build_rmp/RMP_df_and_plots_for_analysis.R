### Build 3 ###

df <- sweep_objs$data_00039_AD0$df %>%
  dplyr::arrange(epoch_t) %>%
  dplyr::mutate(
    state_change = pulse_state != dplyr::lag(pulse_state, default = first(pulse_state)),
    pulse_run    = cumsum(state_change),
    pulse_on_idx  = if_else(pulse_state == "pulse_on",
                            cumsum(state_change & pulse_state == "pulse_on"),
                            NA_integer_),
    pulse_off_idx = if_else(pulse_state == "pulse_off",
                            cumsum(state_change & pulse_state == "pulse_off"),
                            NA_integer_),
    pulse_id = dplyr::case_when(
      pulse_state == "pulse_on"  ~ sprintf("on_%03d",  pulse_on_idx),
      pulse_state == "pulse_off" ~ sprintf("off_%03d", pulse_off_idx)
    )
  ) %>%
  dplyr::group_by(pulse_id) %>%
  dplyr::mutate(
    t0 = min(epoch_t),
    t_rel = epoch_t - t0
  ) %>%
  dplyr::ungroup()

library(dplyr)

seg_info <- df %>%
  dplyr::group_by(pulse_id, pulse_state) %>%
  dplyr::summarise(t_start = min(epoch_t), .groups = "drop") %>%
  dplyr::arrange(t_start)

last_row <- seg_info %>% dplyr::slice_tail(n = 1)

if (last_row$pulse_state == "pulse_on") {
  drop_ids <- seg_info %>% dplyr::slice_tail(n = 2) %>% dplyr::pull(pulse_id)
} else {
  drop_ids <- last_row$pulse_id
}

df <- df %>% dplyr::filter(!pulse_id %in% drop_ids)


pulse_on_summary <- df %>%
  dplyr::filter(pulse_state == "pulse_on") %>%
  dplyr::group_by(pulse_id) %>%
  dplyr::summarise(
    mean_on = mean(tpRMP, na.rm = TRUE),
    peak_neg = min(tpRMP, na.rm = TRUE),
    t_peak   = t_rel[which.min(tpRMP)],
    duration = max(t_rel),
    delta_mean_minus_peak = mean_on - peak_neg,  # (mV) >0 if peak is more negative than mean
    .groups = "drop"
  )

pulse_tau <- df %>%
  dplyr::filter(pulse_state == "pulse_on") %>%
  dplyr::group_by(pulse_id) %>%
  dplyr::summarise(
    peak_tpRMP = min(tpRMP),
    steady_tpRMP = mean(tpRMP[t_rel > 0.8 * max(t_rel)], na.rm = TRUE),
    target = peak_tpRMP + 0.63 * (steady_tpRMP - peak_tpRMP),
    tau = t_rel[which.min(abs(tpRMP - target))],
    .groups = "drop"
  )

pulse_off_summary <- df %>%
  dplyr::filter(pulse_state == "pulse_off") %>%
  dplyr::group_by(pulse_id) %>%
  dplyr::summarise(
    mean_off = mean(tpRMP, na.rm = TRUE),
    sd_off   = sd(tpRMP, na.rm = TRUE),
    duration = max(t_rel),

    # early peak in first 100 ms
    early_peak_off = {
      d100 <- tpRMP[t_rel <= 0.1]
      if (length(d100) == 0) NA_real_ else max(d100, na.rm = TRUE)
    },
    early_t_peak_off = {
      idx <- which(t_rel <= 0.1)
      if (length(idx) == 0) NA_real_ else t_rel[idx][which.max(tpRMP[idx])]
    },

    .groups = "drop"
  )

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

p_all


df %>%
  ggplot(., aes(x = epoch_t, y = tpRMP, color = pulse_state)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  xlim(0,.2) +
  labs(
    x = "Time (s)",
    y = "tpRMP (mV)",
    color = "Pulse state"
  ) +
  theme_classic() +
  theme(legend.position = "none")

df %>%
  ggplot(., aes(x = epoch_t, y = tpRMP, color = pulse_state)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  xlim(0,0.99) +
  labs(
    x = "Time (s)",
    y = "tpRMP (mV)",
    color = "Pulse state"
  ) +
  theme_classic()


