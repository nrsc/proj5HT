## 6) Summary visuals
print(plot_subclass_response_bars(per_cell_plot))
print(plot_mid_late_dot(per_cell_plot, metric = "late_mean_pct"))
print(plot_split_violins(per_cell_plot,
                         subclasses = c("L23_IT","L5_ET","L5_IT"),
                         baseline = 100, gap = 4))

exc_df <- per_cell2 %>%
  dplyr::filter(direction == "excitation") %>%
  dplyr::filter(is.finite(exc_mag), exc_mag > 0) %>%
  dplyr::group_by(blockers_bin) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_exc = mean(exc_mag, na.rm = TRUE),
    sd_exc   = sd(exc_mag, na.rm = TRUE),
    se_exc   = sd_exc / sqrt(n),
    ci95     = 1.96 * se_exc,
    .groups = "drop"
  )

p_exc_bar <- ggplot(exc_df, aes(x = blockers_bin, y = mean_exc)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_exc - ci95, ymax = mean_exc + ci95), width = 0.2, linewidth = 0.6) +
  geom_text(aes(label = paste0("n=", n)), vjust = -0.6, size = 3.6) +
  labs(
    x = NULL,
    y = "Excitation magnitude (% above baseline)",
    title = "Excitation magnitude vs blockers",
    subtitle = "Bars show mean ± 95% CI (excitatory cells only)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_exc_bar)

p_exc_dist <- per_cell2 %>%
  dplyr::filter(direction == "excitation", is.finite(exc_mag), exc_mag > 0) %>%
  ggplot(aes(x = blockers_bin, y = exc_mag)) +
  geom_violin(trim = TRUE, alpha = 0.5) +
  geom_jitter(width = 0.12, alpha = 0.35, size = 1) +
  labs(
    x = NULL,
    y = "Excitation magnitude (% above baseline)",
    title = "Excitation magnitude distribution (QC)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

print(p_exc_dist)


## 7) Diagnose “insufficient data”
diag <- diagnose_insufficient(df_prepped, baseline_min_pts = 5, post_min_pts = 20, post_window = c(0,60))
print(table(diag$problem))
