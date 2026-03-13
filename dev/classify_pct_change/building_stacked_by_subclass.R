suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(tidyr)
})

plot_biphasic_stacked_by_subclass <- function(df_prepped,
                                              per_cell,
                                              biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                              xlim = c(-5, 20),
                                              k_smooth = 11,
                                              stack_gap = 35,
                                              rapid_window = c(0, 10),
                                              show_rapid_points = TRUE,
                                              title = NULL) {

  biphasic_order_keep <- match.arg(biphasic_order_keep)

  # --- choose biphasic cells to keep (ONLY bring cell_name to avoid join suffixes) ---
  keep_cells <- per_cell %>%
    dplyr::filter(
      response_type == "biphasic",
      biphasic_order == biphasic_order_keep
    ) %>%
    dplyr::pull(cell_name) %>%
    unique()

  if (length(keep_cells) == 0) stop("No biphasic cells found for ", biphasic_order_keep)

  d_stack <- df_prepped %>%
    dplyr::filter(cell_name %in% keep_cells) %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2])

  # --- robust subclass column selection after any prior joins ---
  # Prefer plain 'assigned_subclass' if present; otherwise look for suffix variants.
  subclass_candidates <- c("assigned_subclass", "assigned_subclass.x", "assigned_subclass.y")
  subclass_found <- subclass_candidates[subclass_candidates %in% names(d_stack)][1]
  if (is.na(subclass_found)) {
    stop("df_prepped does not contain assigned_subclass (or .x/.y). Add it before plotting.")
  }

  d_stack <- d_stack %>%
    dplyr::mutate(subclass_facet = .data[[subclass_found]])

  # --- smoothing per trace within subclass ---
  d_stack <- d_stack %>%
    dplyr::group_by(subclass_facet, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  # stacking order inside each subclass facet
  ord <- d_stack %>%
    dplyr::distinct(subclass_facet, cell_name) %>%
    dplyr::arrange(subclass_facet, cell_name) %>%
    dplyr::group_by(subclass_facet) %>%
    dplyr::mutate(trace_idx = dplyr::row_number()) %>%
    dplyr::ungroup()

  d_stack <- d_stack %>%
    dplyr::left_join(ord, by = c("subclass_facet","cell_name")) %>%
    dplyr::mutate(
      offset   = (trace_idx - 1) * stack_gap,
      y_raw    = pct_of_baseline + offset,
      y_smooth = pct_smooth + offset
    )

  base_lines <- d_stack %>%
    dplyr::distinct(subclass_facet, cell_name, trace_idx) %>%
    dplyr::mutate(
      offset = (trace_idx - 1) * stack_gap,
      y0 = 100 + offset
    )

  # --- rapid markers (optional) ---
  marks <- NULL
  if (isTRUE(show_rapid_points)) {
    marks <- per_cell %>%
      dplyr::filter(
        response_type == "biphasic",
        biphasic_order == biphasic_order_keep
      ) %>%
      dplyr::select(cell_name, rapid_min_time, rapid_min_pct, rapid_max_time, rapid_max_pct) %>%
      tidyr::pivot_longer(
        cols = dplyr::starts_with("rapid_"),
        names_to = c("which","field"),
        names_pattern = "rapid_(min|max)_(time|pct)",
        values_to = "value"
      ) %>%
      tidyr::pivot_wider(names_from = field, values_from = value) %>%
      dplyr::rename(mark = which) %>%
      dplyr::filter(is.finite(time), is.finite(pct)) %>%
      dplyr::left_join(
        d_stack %>% dplyr::distinct(cell_name, subclass_facet, trace_idx),
        by = "cell_name"
      ) %>%
      dplyr::mutate(
        offset = (trace_idx - 1) * stack_gap,
        y = pct + offset,
        mark = factor(mark, levels = c("max","min"))
      )
  }

  # Colors: RED = excitation, BLUE = inhibition
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  if (is.null(title)) {
    title <- if (biphasic_order_keep == "exc_then_inh")
      "Rapid biphasic (excitation → inhibition), stacked by subclass"
    else
      "Rapid biphasic (inhibition → excitation), stacked by subclass"
  }

  ggplot(d_stack, aes(x = time)) +
    facet_grid(subclass_facet ~ ., scales = "free_y") +
    geom_hline(data = base_lines, aes(yintercept = y0),
               linetype = 2, linewidth = 0.35, color = "grey40") +
    annotate("rect",
             xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.06, fill = "grey30") +
    geom_line(aes(y = y_raw, group = cell_name),
              color = "grey70", linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = y_smooth, group = cell_name),
              color = "black", linewidth = 0.6, na.rm = TRUE) +
    { if (!is.null(marks)) geom_point(
      data = marks,
      aes(x = time, y = y, shape = mark, color = mark),
      size = 2.2
    ) else NULL } +
    scale_shape_manual(values = c(max = 24, min = 25)) +
    scale_color_manual(values = c(max = col_exc, min = col_inh)) +
    labs(
      title = title,
      subtitle = "Subclass taken from df_prepped; baseline=100; rapid markers shown as triangles",
      x = "Time (s)",
      y = "% baseline + offset",
      color = "rapid feature",
      shape = "rapid feature"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
}

p_excfirst_sub <- plot_biphasic_stacked_by_subclass(
  df_prepped = df_prepped,
  per_cell   = per_cell2,
  biphasic_order_keep = "exc_then_inh",
  xlim = c(-5, 20),
  k_smooth = 2,
  stack_gap = 0
)
print(p_excfirst_sub)

