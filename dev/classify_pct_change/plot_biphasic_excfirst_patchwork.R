suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------
# Patchwork grid of biphasic (rapid) responses: exc -> inh
# Requires:
# - per_cell2 has: cell_name, response_type, biphasic_order, assigned_subclass
# - df_prepped is your prepped trace table
# - plot_bucket_qc(df_prepped, cell_id, ...) exists (from your pipeline)
# ---------------------------------------------------------

plot_biphasic_excfirst_patchwork <- function(df_prepped,
                                             per_cell2,
                                             subclasses = NULL,     # e.g. c("L23_IT","L5_ET","L5_IT")
                                             n_max = 24,            # max panels to show
                                             ncol = 4,              # grid columns
                                             k_smooth = 11,
                                             smooth_method = "rollmean",
                                             min_excursion_pct = 5,
                                             rapid_window = c(0,10),
                                             mid_window = c(10,40),
                                             late_slope_window = c(36,56),
                                             late_mean_window  = c(46,56),
                                             post_window = c(0,60)) {

  stopifnot("cell_name" %in% names(per_cell2))
  stopifnot("response_type" %in% names(per_cell2))
  stopifnot("biphasic_order" %in% names(per_cell2))

  pick <- per_cell2 %>%
    dplyr::filter(response_type == "biphasic",
                  biphasic_order == "exc_then_inh")

  if (!is.null(subclasses)) {
    if (!"assigned_subclass" %in% names(pick)) stop("assigned_subclass missing in per_cell2")
    pick <- pick %>% dplyr::filter(assigned_subclass %in% subclasses)
  }

  if (nrow(pick) == 0) stop("No biphasic exc_then_inh cells found with the current filters.")

  # pick up to n_max cells (optionally bias toward strongest peaks if available)
  if ("response_peak" %in% names(pick)) {
    pick <- pick %>% dplyr::arrange(dplyr::desc(response_peak))
  }

  pick <- pick %>% dplyr::slice_head(n = n_max)

  plots <- lapply(pick$cell_name, function(cid) {
    plot_bucket_qc(
      df_prepped,
      cell_id = cid,
      k_smooth = k_smooth,
      smooth_method = smooth_method,
      min_excursion_pct = min_excursion_pct,
      rapid_window = rapid_window,
      mid_window = mid_window,
      late_slope_window = late_slope_window,
      late_mean_window = late_mean_window,
      post_window = post_window
    ) +
      theme(
        plot.title = element_text(size = 9),
        plot.subtitle = element_text(size = 7),
        axis.title = element_blank()
      )
  })

  wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = paste0("Rapid biphasic responses (excitation → inhibition) | n = ", nrow(pick)),
      theme = theme(plot.title = element_text(face = "bold"))
    )
}

# ---------------------------
# Example usage
# ---------------------------

p_patch <- plot_biphasic_excfirst_patchwork(
  df_prepped = df_prepped,
  per_cell2  = per_cell2,
  subclasses = c("L23_IT","L5_ET","L5_IT"),  # optional
  n_max = 20,
  ncol = 4,
  k_smooth = 11,
  smooth_method = "rollmean",
  min_excursion_pct = 5
)

print(p_patch)

# Save for slides:
#ggsave("figs/qc/biphasic_exc_first_patchwork.png", p_patch, width = 14, height = 10, dpi = 250)
