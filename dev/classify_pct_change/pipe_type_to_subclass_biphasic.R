suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(tidyr)
  library(rlang)
})

make_blockers_bin <- function(x) {
  dplyr::if_else(is.na(x) | x == "", "no_blockers", "has_blockers")
}

prep_biphasic_stack_df_facet <- function(df_prepped,
                                         per_cell,
                                         biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                         facet_var = "assigned_subclass",
                                         xlim = c(-5, 20),
                                         k_smooth = 11) {

  # ---- safety checks ----
  need_df <- c("cell_name","time","pct_of_baseline")
  miss_df <- setdiff(need_df, names(df_prepped))
  if (length(miss_df) > 0) stop("df_prepped missing: ", paste(miss_df, collapse=", "))

  need_pc <- c("cell_name","response_type","biphasic_order",
               "rapid_min_time","rapid_min_pct","rapid_max_time","rapid_max_pct")
  miss_pc <- setdiff(need_pc, names(per_cell))
  if (length(miss_pc) > 0) stop("per_cell missing: ", paste(miss_pc, collapse=", "))

  if (!facet_var %in% names(df_prepped)) {
    stop("facet_var '", facet_var, "' not found in df_prepped. Available: ",
         paste(names(df_prepped), collapse = ", "))
  }

  # keep only rapid biphasic cells of requested order
  keep_cells <- per_cell %>%
    dplyr::filter(.data[["response_type"]] == "biphasic",
                  .data[["biphasic_order"]] == biphasic_order_keep) %>%
    dplyr::pull(.data[["cell_name"]]) %>%
    unique()

  if (length(keep_cells) == 0) {
    stop("No biphasic cells found with biphasic_order == ", biphasic_order_keep)
  }

  d <- df_prepped %>%
    dplyr::filter(.data[["cell_name"]] %in% keep_cells) %>%
    dplyr::mutate(
      time = as.numeric(.data[["time"]]),
      blockers_bin = if ("blockers" %in% names(df_prepped)) make_blockers_bin(.data[["blockers"]]) else NA_character_
    ) %>%
    dplyr::filter(.data[["time"]] >= xlim[1], .data[["time"]] <= xlim[2]) %>%
    dplyr::group_by(.data[[facet_var]], .data[["cell_name"]]) %>%
    dplyr::arrange(.data[["time"]], .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(.data[["pct_of_baseline"]],
                                 k = k_smooth, fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  # stacking order within each facet
  ord <- d %>%
    dplyr::distinct(.data[[facet_var]], .data[["cell_name"]]) %>%
    dplyr::arrange(.data[[facet_var]], .data[["cell_name"]]) %>%
    dplyr::group_by(.data[[facet_var]]) %>%
    dplyr::mutate(trace_idx = dplyr::row_number()) %>%
    dplyr::ungroup()

  dplyr::left_join(d, ord, by = c(facet_var, "cell_name"))
}

prep_rapid_markers_facet <- function(d_stack,
                                     per_cell,
                                     facet_var,
                                     biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                     stack_gap = 35) {

  off <- d_stack %>%
    dplyr::distinct(.data[[facet_var]], .data[["cell_name"]], .data[["trace_idx"]]) %>%
    dplyr::mutate(offset = (.data[["trace_idx"]] - 1) * stack_gap)

  marks <- per_cell %>%
    dplyr::filter(.data[["response_type"]] == "biphasic",
                  .data[["biphasic_order"]] == biphasic_order_keep) %>%
    dplyr::select(.data[["cell_name"]],
                  .data[["rapid_min_time"]], .data[["rapid_min_pct"]],
                  .data[["rapid_max_time"]], .data[["rapid_max_pct"]]) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("rapid_"),
      names_to = c("which","field"),
      names_pattern = "rapid_(min|max)_(time|pct)",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(names_from = field, values_from = value) %>%
    dplyr::rename(mark = which) %>%
    dplyr::filter(is.finite(.data[["time"]]), is.finite(.data[["pct"]])) %>%
    dplyr::left_join(off, by = "cell_name") %>%
    dplyr::mutate(
      y = .data[["pct"]] + .data[["offset"]],
      mark = factor(.data[["mark"]], levels = c("max","min"))
    )

  marks
}

prep_baseline_lines_facet <- function(d_stack, facet_var, stack_gap = 35) {
  d_stack %>%
    dplyr::distinct(.data[[facet_var]], .data[["cell_name"]], .data[["trace_idx"]]) %>%
    dplyr::mutate(
      offset = (.data[["trace_idx"]] - 1) * stack_gap,
      y0 = 100 + offset
    )
}

plot_biphasic_stacked_facet <- function(df_prepped,
                                        per_cell,
                                        biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                        facet_var = "assigned_subclass",
                                        xlim = c(-5, 20),
                                        k_smooth = 11,
                                        stack_gap = 35,
                                        rapid_window = c(0, 10),
                                        show_rapid_points = TRUE,
                                        title = NULL) {

  d_stack <- prep_biphasic_stack_df_facet(
    df_prepped, per_cell,
    biphasic_order_keep = biphasic_order_keep,
    facet_var = facet_var,
    xlim = xlim,
    k_smooth = k_smooth
  )

  d_stack <- d_stack %>%
    dplyr::mutate(
      offset   = (.data[["trace_idx"]] - 1) * stack_gap,
      y_raw    = .data[["pct_of_baseline"]] + .data[["offset"]],
      y_smooth = .data[["pct_smooth"]] + .data[["offset"]]
    )

  base_lines <- prep_baseline_lines_facet(d_stack, facet_var = facet_var, stack_gap = stack_gap)

  marks <- NULL
  if (isTRUE(show_rapid_points)) {
    marks <- prep_rapid_markers_facet(
      d_stack, per_cell,
      facet_var = facet_var,
      biphasic_order_keep = biphasic_order_keep,
      stack_gap = stack_gap
    )
  }

  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  if (is.null(title)) {
    title <- if (biphasic_order_keep == "exc_then_inh") {
      "Rapid biphasic responses (excitation → inhibition)"
    } else {
      "Rapid biphasic responses (inhibition → excitation)"
    }
  }

  ggplot(d_stack, aes(x = .data[["time"]])) +
    facet_grid(rows = vars(!!sym(facet_var)), scales = "free_y") +

    geom_hline(data = base_lines, aes(yintercept = .data[["y0"]]),
               linetype = 2, linewidth = 0.35, color = "grey40") +

    annotate("rect",
             xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.06, fill = "grey30") +

    geom_line(aes(y = .data[["y_raw"]], group = .data[["cell_name"]]),
              color = "grey70", linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = .data[["y_smooth"]], group = .data[["cell_name"]]),
              color = "black", linewidth = 0.6, na.rm = TRUE) +

    { if (!is.null(marks)) geom_point(
      data = marks,
      aes(x = .data[["time"]], y = .data[["y"]], shape = .data[["mark"]], color = .data[["mark"]]),
      size = 2.2
    ) else NULL } +

    scale_shape_manual(values = c(max = 24, min = 25)) +
    scale_color_manual(values = c(max = col_exc, min = col_inh)) +

    labs(
      title = title,
      subtitle = paste0("Baseline = mean[-5,0] → 100%; facet = ", facet_var),
      x = "Time (s)",
      y = "% baseline (100) + vertical offset",
      shape = "rapid feature",
      color = "rapid feature"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey70")
    )
}

p_excfirst_subclass <- plot_biphasic_stacked_facet(
  df_prepped = df_prepped,
  per_cell   = per_cell2,
  biphasic_order_keep = "exc_then_inh",
  facet_var = "assigned_subclass",
  xlim = c(-5, 30),
  k_smooth = 5,
  stack_gap = 0,
  rapid_window = c(0, 10),
  show_rapid_points = TRUE
)

print(p_excfirst_subclass)

names(per_cell2)
table(per_cell2$response_type, useNA="ifany")
table(per_cell2$biphasic_order, useNA="ifany")

