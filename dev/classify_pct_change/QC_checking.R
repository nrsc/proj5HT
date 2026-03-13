
df_cv <- df %>%
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
  dplyr::filter(assigned_subclass %in% c("L23_IT","L5_ET","L5_IT")) %>%
  # ---- optional: restrict window used for CV/rate (recommended) ----
  dplyr::filter(time >= -20, time <= 0) %>%
# ---- compute per-cell spiking variability features ----
dplyr::group_by(cell_name) %>%
  dplyr::mutate(
    instRate_mean = mean(instRate, na.rm = TRUE),
    instRate_sd   = sd(instRate, na.rm = TRUE),
    instRate_cv   = dplyr::if_else(instRate_mean > 0,
                                   instRate_sd / instRate_mean,
                                   NA_real_),
    instRate_n    = sum(!is.na(instRate))
  ) %>%
  dplyr::ungroup()

df_cv <- df_cv %>%
  dplyr::mutate(
    bad_spiking = instRate_mean < 2 & instRate_cv > 0.5
  )

blCV <- df_cv %>%
  dplyr::distinct(cell_name, assigned_subclass, instRate_mean, instRate_cv)
fullCV <- df_use %>%
  dplyr::distinct(cell_name, assigned_subclass, instRate_mean, instRate_cv)

cv_cells[which(cv_cells$instRate_cv > 0.6),]

p1 = df_cv %>%  dplyr::distinct(cell_name, assigned_subclass, instRate_mean, instRate_cv) %>%
ggplot(., aes(x = instRate_mean, y = instRate_cv, colour = assigned_subclass)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  theme_bw() +
  labs(title = "CV vs instRate for baseline -20:0", x = "Mean instRate (Hz, log10)", y = "CV(instRate)", colour = "Subclass")

p2 = df_use %>%  dplyr::distinct(cell_name, assigned_subclass, instRate_mean, instRate_cv) %>%
ggplot(., aes(x = instRate_mean, y = instRate_cv, colour = assigned_subclass)) +
  geom_point(alpha = 0.8) +
  scale_x_log10() +
  theme_bw() +
  labs(title = "CV vs instRate for all isi", x = "Mean instRate (Hz, log10)", y = "CV(instRate)", colour = "Subclass")

grid.arrange(p1,p2, ncol = 2)





### Low spike rate low cv
srt = load5HT("QN26.26.004.20.01.01", rehydrate = FALSE)
plotBasicPuff(srt, xlim = c(-40,50))
out <- plot_sweeps_py(
  nwb_path = srt$files$nwb,
  sweeps = c(
    srt$dfs$spikeTTL$sweep_names[["baseline"]],
    srt$dfs$spikeTTL$sweep_names[["spikeTTL"]]
  ),
  fmt = "pdf",
  title = paste(srt$cell, unique(srt$dfs$spikeTTL$spike_puff_output$assigned_subclass), sep = " -- "),
  re_zero = TRUE,
  gap_s = 0,

  ttl_time = srt$dfs$spikeTTL$ttlTime,
  ttl_on_sweep_index = 1,

  pub_void_axes = TRUE,
  y_unit_target = "mV",
  scalebar_x = 10,
  scalebar_label_x = "10 s",
  scalebar_y = 20,
  scalebar_label_y = "20 mV"
)

message("Saved: ", out)

rstudioapi::viewer(out)


## High CV, Mid frequency
srt = load5HT("QN25.26.014.11.18.01", rehydrate = FALSE)
plotBasicPuff(srt, xlim = c(-30,40))
out <- plot_sweeps_py(
  nwb_path = srt$files$nwb,
  sweeps = c(
    srt$dfs$spikeTTL$sweep_names[["baseline"]],
    srt$dfs$spikeTTL$sweep_names[["spikeTTL"]]
  ),
  fmt = "pdf",

  re_zero = TRUE,
  gap_s = 0,

  ttl_time = srt$dfs$spikeTTL$ttlTime,
  ttl_on_sweep_index = 1,

  pub_void_axes = TRUE,
  y_unit_target = "mV",
  scalebar_x = 10,
  scalebar_label_x = "10 s",
  scalebar_y = 20,
  scalebar_label_y = "20 mV"
)

message("Saved: ", out)

rstudioapi::viewer(out)




