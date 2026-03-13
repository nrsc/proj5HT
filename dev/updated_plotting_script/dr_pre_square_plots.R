df = data.frame(comp5HT$srtPuff)
df$expCon = gsub("standard_puff", "Standard_Puff", df$expCon)
df$assigned_subclass = gsub("/", "", df$assigned_subclass)
df$Date = as.Date(df$Date)

df = df %>% mutate(
  x_bin = floor(time / 5) * 5 + 5
)
df = df %>% dplyr::filter(
  !cell_name %in% c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03"),
)

df2 = df %>%
  dplyr::filter(
    Date >= as.Date("2025-09-17") &
      Date <= as.Date("2026-01-14")
  )
df0 = df %>%
  dplyr::filter(
    !dplyr::between(
      Date,
      as.Date("2025-09-17"),
      as.Date("2026-01-14")
    ))

df_sel = df0 %>% dplyr::filter(
  bath == "none",
  assigned_subclass %in% c("L23_IT", "L5_ET", "L5_IT"),
  grepl("^5HT", puff)
)
#

df_prepped <- prep_puff_df(df_sel, baseline_window = c(-5, 0))

per_cell <- extract_bucket_features(df_prepped, min_excursion_pct = 10, k_smooth = 11, smooth_method = "rollmean")
per_cell2 <- per_cell %>% add_biphasic_features(min_excursion_pct = 10) # implement below or reuse your existing helper
per_cell2 <- per_cell2 %>%
  add_response_call(mild = 10, moderate = 25, strong = 50)

plot_bucket_qc(df_prepped, "QN26.26.004.20.03.03", k_smooth = 11, smooth_method = "rollmean")
# print(p_qc)
#
# note: I did not automatically re-implement your 'add_biphasic_features' here verbatim;
# you can either reuse the original function you provided or call the short one below.

# ---- panel figures ---------------
plot_biphasic_modalities_by_subclass(df_prepped, per_cell2,
                                             xlim = c(-5, 30),
                                             mean_smooth = "spline",
                                             spline_df = 50, k_smooth = 11, alpha_traces = 0.3, trace_color = "purple1")



per_cell2 %>% dplyr::filter(assigned_subclass == "L5_IT", direction == "excitation", rapid_pattern == "rapid_excitation_only")

unique(per_cell2$rapid_pattern)

plot_exc_inh_by_subclass(df_prepped, per_cell2,
                                 xlim = c(-5, 30), k_smooth = 11,
                                 min_cells_per_time = 3,
                                 mean_smooth = "spline",
                                 spline_df = 25,
                                 alpha_traces = .3,
                                 trace_color = "purple1",
                                 mean_color = "black",
                                 y_scale = 20
                                 )



# END OF SCRIPT


plot_df = prep_plot_df(per_cell2)
plot_centered_violin(plot_df)
plot_stacked_100(plot_df)
plot_heatmap_prop(per_cell2)

