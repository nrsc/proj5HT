# source("puff_helpers.R")
# source("puff_processing.R")
# source("puff_plotting.R")

# df is your original long dataframe
# names(df) already include x_bin with explicit zeros for no-spike bins

df_prepped <- prepare_puff_df(
  df,
  baseline_window = c(-5, 0),
  subclass_keep = c("L23_IT", "L5_ET", "L5_IT"),
  bath_keep = "none",
  puff_pattern = "^5HT"
)

per_cell <- build_per_cell_classifier(
  df_prepped,
  rapid_window = c(0, 10),
  mid_window = c(10, 40),
  late_window = c(40, 60),
  post_window = c(0, 60),
  min_excursion_pct = 10,
  min_persist_s = 0.5,
  smooth_method = "rollmean",
  smooth_k = 3
)

# QC
p1 <- plot_cell_qc(df_prepped, cell_id = per_cell$cell_name[1], smooth_k = 3)
print(p1)

# Summary figures
print(plot_response_composition(per_cell))
print(plot_response_heatmap(per_cell))
print(plot_centered_violin(per_cell))
print(plot_exc_inh_traces(df_prepped, per_cell))
