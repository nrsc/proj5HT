source("~/proj5HT/claude/checking_consistency.r")
## 1) Prep
#df_prepped <- prep_puff_df(df0, baseline_window = c(-5, 0))

#df_prepped <- prep_puff_df(df0)

df_prepped <- prep_puff_df(df0, add_xbin_anchors = TRUE, x_bin_s = 2, xlim = c(-10, 50),
                    anchor_id_cols = c("cell_name"))

# plot_bucket_qc(df_prepped, "QN25.26.014.11.12.01",
#                k_smooth = 7,
#                smooth_method = "rollmean",
#                rapid_window = c(0, 10),
#                mid_window = c(10, 30),
#                late_slope_window = c(30, 50),  # last 20s
#                late_mean_window  = c(45, 50),  # last 10s
#                post_window = c(0, 60),
# )

## 2) Extract per-cell bucket features
per_cell <- extract_bucket_features(df_prepped,
                                    rapid_window = c(0, 15),
                                    mid_window = c(10, 35),
                                    late_slope_window = c(40, 50),
                                    late_mean_window  = c(40, 50),
                                    post_window = c(0, 35),
                                    min_excursion_pct = 15,
                                    k_smooth = 5,
                                    min_persist_s = 0.5,
                                    smooth_method = "rollmean")

per_cell = extract_bucket_features(
  df_prepped,
  axis = "x_bin",
  rapid_window = c(0, 10),
  post_window  = c(0, 60),
  k_smooth_bins = 3
)


# per_cell2 <- per_cell %>%
#   add_rapid_zero_inh(
#     df_prepped,
#     rapid_window = c(0, 10),
#     rate_col = "instRate",
#     zero_thr_hz = 0.1,
#     min_zero_s = 3
#   ) %>%
#   add_biphasic_features_simple(
#     min_excursion_pct = 15,
#     min_sep_s = 0.25,
#     use_mid_rebound_gate = TRUE,
#     rebound_mid_thr = 5,
#     still_inh_mid_thr = 10
#   )


per_cell2 <- per_cell %>%
  add_biphasic_features_simple(
    min_excursion_pct = 15,
    min_sep_s = 0.25,
    use_mid_rebound_gate = TRUE,
    rebound_mid_thr = 5,
    still_inh_mid_thr = 10
  ) %>%
  dplyr::mutate(
    response_type  = factor(response_type, levels = c("excitation","inhibition","biphasic","no_change")),
    biphasic_order = factor(biphasic_order, levels = c("exc_then_inh","inh_then_exc","ambiguous"))
  )


per_cell2 <- per_cell %>%
  add_biphasic_features_consistent(
    min_excursion_pct = 15,
    rescue_biphasic   = TRUE,
    rescue_inh_window = c(-1, 10),
    rescue_area_thr   = 1,
    min_sep_s         = 2,
    rel_inh_to_exc    = 0.35,
    deep_dip_pct      = 70,
    rebound_mid_thr       = 5,
    raw_rebound_mid_thr   = 2,
    raw_use = "mean",   # try "max" if you want to be more permissive
    still_inh_mid_thr     = 10
  )  %>%
  dplyr::mutate(
    response_type  = as.character(response_type),
    biphasic_order = as.character(biphasic_order)
  )


# Checks for improperly labeled excitation
per_cell2 %>%
  dplyr::filter(response_type == "excitation") %>%
  dplyr::mutate(
    has_any_inh = is.finite(abs_min_pct) & abs_min_pct <= (100 - 10)
  ) %>%
  dplyr::count(has_any_inh)

per_cell2 %>%
  dplyr::filter(response_type == "excitation",
                is.finite(abs_min_pct),
                abs_min_pct <= 90) %>%
  dplyr::select(cell_name, abs_min_pct, abs_min_time,
               abs_max_pct, abs_max_time, biphasic_order)


# use:
print(plot_biphasic_first_by_subclass(per_cell2))
print(plot_basic_by_subclass(per_cell2))


p_ei <- plot_exc_inh_by_subclass(
  df_prepped, per_cell2,
  xlim = c(-5, 39),
  k_smooth = 11,
  min_cells_per_time = 2,
  mean_smooth = "spline",
  spline_df = 10,
  alpha_traces = .3,
  trace_color = "purple1",
  mean_color = "black",
  y_scale = 20
)
print(p_ei)

p_bi <-plot_biphasic_modalities_by_subclass(
  df_prepped = df_prepped,
  per_cell = per_cell2,
  xlim = c(-5, 39),
  mean_smooth = "spline",
  spline_df = 20,
  k_smooth = 11,
  alpha_traces = 0.3,
  trace_color = "purple1",
  facet_scales = "fixed",     # <-- MATCH Y SCALES
  y_bal_log = TRUE,
  y_scale = 8                 # keep your log-ish scaling
)

print(p_bi)

p_4 <- plot_4row_modalities_by_subclass(
  df_prepped, per_cell2,
  trace_alpha = 0.4,
  trace_color = "purple1",

  mean_color = "navyblue",
  mean_lwd = 1,

  mean_smooth = "spline",
  spline_df = 50,

  facet_scales = "fixed",
  y_bal_log = TRUE,
  y_scale = 8
)

print(p_4)


## 5) Save a bunch for visual checking
#pick some cells from your target subclasses
plot_bucket_qc_many(df_prepped, per_cell_plot$cell_name,
                    out_dir = "qc_bucket_plots",
                    k_smooth = 3,
                    smooth_method = "rollmean",
                    min_excursion_pct = 5)

# unique(df_prepped$puff)
#
plot_bucket_qc(df_prepped, "QN26.26.004.20.03.02",
               k_smooth = ,
               smooth_method = "rollmean",
               rapid_window = c(0, 10),
               mid_window = c(10, 30),
               late_slope_window = c(30, 50),  # last 20s
               late_mean_window  = c(45, 50),  # last 10s
               post_window = c(0, 60),
)


if (!"blockers" %in% names(per_cell2)) {
  blockers_map <- df %>%
    dplyr::distinct(cell_name, blockers)

  per_cell2 <- per_cell2 %>%
    dplyr::left_join(blockers_map, by = "cell_name")
}

per_cell2 <- per_cell2 %>%
  dplyr::mutate(
    blockers_bin = dplyr::if_else(is.na(blockers) | blockers == "", "no_blockers", "has_blockers"),
    blockers_bin = factor(blockers_bin, levels = c("no_blockers", "has_blockers"))
  )
