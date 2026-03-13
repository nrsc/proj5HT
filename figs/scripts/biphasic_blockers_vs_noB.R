suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(tidyr)
})

# --- helpers -------------------------------------------------

make_blockers_bin <- function(x) {
  dplyr::if_else(is.na(x) | x == "", "no_blockers", "has_blockers")
}

# Build plotting table for stacked traces
# per_cell MUST include:
# cell_name, response_type, biphasic_order,
# rapid_min_time, rapid_min_pct, rapid_max_time, rapid_max_pct
prep_biphasic_stack_df <- function(df_prepped,
                                   per_cell,
                                   biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                   xlim = c(-5, 20),
                                   k_smooth = 11) {

  biphasic_order_keep <- match.arg(biphasic_order_keep)

  stopifnot(all(c("cell_name","time","pct_of_baseline","blockers") %in% names(df_prepped)))
  stopifnot(all(c("cell_name","response_type","biphasic_order",
                  "rapid_min_time","rapid_min_pct",
                  "rapid_max_time","rapid_max_pct") %in% names(per_cell)))

  keep_cells <- per_cell %>%
    dplyr::filter(
      response_type == "biphasic",
      biphasic_order == biphasic_order_keep
    ) %>%
    dplyr::pull(cell_name) %>%
    unique()

  if (length(keep_cells) == 0) {
    stop("No biphasic cells found with biphasic_order == ", biphasic_order_keep)
  }

  d <- df_prepped %>%
    dplyr::filter(cell_name %in% keep_cells) %>%
    dplyr::mutate(
      blockers_bin = make_blockers_bin(blockers),
      time = as.numeric(time)
    ) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(blockers_bin, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth,
                                 fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  # stacking order inside each facet
  ord <- d %>%
    dplyr::distinct(blockers_bin, cell_name) %>%
    dplyr::arrange(blockers_bin, cell_name) %>%
    dplyr::group_by(blockers_bin) %>%
    dplyr::mutate(trace_idx = dplyr::row_number()) %>%
    dplyr::ungroup()

  d %>% dplyr::left_join(ord, by = c("blockers_bin","cell_name"))
}

# Rapid markers table (min/max in rapid window), aligned to stacked offsets
prep_rapid_markers <- function(d_stack,
                               per_cell,
                               biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                               stack_gap = 35) {
  biphasic_order_keep <- match.arg(biphasic_order_keep)

  off <- d_stack %>%
    dplyr::distinct(blockers_bin, cell_name, trace_idx) %>%
    dplyr::mutate(offset = (trace_idx - 1) * stack_gap)

  marks <- per_cell %>%
    dplyr::filter(
      response_type == "biphasic",
      biphasic_order == biphasic_order_keep
    ) %>%
    dplyr::select(cell_name,
                  rapid_min_time, rapid_min_pct,
                  rapid_max_time, rapid_max_pct) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("rapid_"),
      names_to = c("which","field"),
      names_pattern = "rapid_(min|max)_(time|pct)",
      values_to = "value"
    ) %>%
    tidyr::pivot_wider(names_from = field, values_from = value) %>%
    dplyr::rename(mark = which) %>%
    dplyr::filter(is.finite(time), is.finite(pct)) %>%
    dplyr::left_join(off, by = "cell_name") %>%
    dplyr::mutate(
      y = pct + offset,
      mark = factor(mark, levels = c("max","min"))
    )

  marks
}

prep_baseline_lines <- function(d_stack, stack_gap = 35) {
  d_stack %>%
    dplyr::distinct(blockers_bin, cell_name, trace_idx) %>%
    dplyr::mutate(
      offset = (trace_idx - 1) * stack_gap,
      y0 = 100 + offset
    )
}

# --- main stacked plot ---------------------------------------

plot_biphasic_stacked <- function(df_prepped,
                                  per_cell,
                                  biphasic_order_keep = c("exc_then_inh","inh_then_exc"),
                                  xlim = c(-5, 20),
                                  k_smooth = 11,
                                  stack_gap = 35,
                                  rapid_window = c(0, 10),
                                  show_rapid_points = TRUE,
                                  title = NULL) {

  biphasic_order_keep <- match.arg(biphasic_order_keep)

  d_stack <- prep_biphasic_stack_df(
    df_prepped, per_cell,
    biphasic_order_keep = biphasic_order_keep,
    xlim = xlim,
    k_smooth = k_smooth
  )

  d_stack <- d_stack %>%
    dplyr::mutate(
      offset   = (trace_idx - 1) * stack_gap,
      y_raw    = pct_of_baseline + offset,
      y_smooth = pct_smooth + offset
    )

  base_lines <- prep_baseline_lines(d_stack, stack_gap = stack_gap)

  marks <- NULL
  if (isTRUE(show_rapid_points)) {
    marks <- prep_rapid_markers(
      d_stack, per_cell,
      biphasic_order_keep = biphasic_order_keep,
      stack_gap = stack_gap
    )
  }

  # RED = excitation, BLUE = inhibition
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  if (is.null(title)) {
    title <- if (biphasic_order_keep == "exc_then_inh") {
      "Rapid biphasic responses (excitation → inhibition)"
    } else {
      "Rapid biphasic responses (inhibition → excitation)"
    }
  }

  ggplot(d_stack, aes(x = time)) +
    facet_grid(blockers_bin ~ ., scales = "free_y") +

    # per-trace baseline at 100+offset
    geom_hline(
      data = base_lines, aes(yintercept = y0),
      linetype = 2, linewidth = 0.35, color = "grey40"
    ) +

    # rapid window shading
    annotate("rect",
             xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.06, fill = "grey30") +

    # raw + smoothed traces
    geom_line(aes(y = y_raw, group = cell_name),
              color = "grey70", linewidth = 0.25, alpha = 0.6) +
    geom_line(aes(y = y_smooth, group = cell_name),
              color = "black", linewidth = 0.6, na.rm = TRUE) +

    # rapid biphasic points only
    { if (!is.null(marks)) geom_point(
      data = marks,
      aes(x = time, y = y, shape = mark, color = mark),
      size = 2.2
    ) else NULL } +

    scale_shape_manual(values = c(max = 24, min = 25),
                       labels = c(max = "rapid_max", min = "rapid_min")) +
    scale_color_manual(values = c(max = col_exc, min = col_inh),
                       labels = c(max = "rapid_max", min = "rapid_min")) +

    labs(
      title = title,
      subtitle = "Baseline = mean[-5,0] → 100%; stacked by blockers (has vs none)",
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

# --- usage ----------------------------------------------------
# Excitation-first stacked:
p_excfirst <- plot_biphasic_stacked(
  df_prepped = df_prepped,
  per_cell   = per_cell2,
  biphasic_order_keep = "exc_then_inh",
  xlim = c(-5, 30),
  k_smooth = 5,
  stack_gap = 0,
  rapid_window = c(0, 10),
  show_rapid_points = TRUE
)

# Inhibition-first stacked (NEW):
p_inhfirst <- plot_biphasic_stacked(
  df_prepped = df_prepped,
  per_cell   = per_cell2,
  biphasic_order_keep = "inh_then_exc",
  xlim = c(-5, 30),
  k_smooth = 5,
  stack_gap = 0,
  rapid_window = c(0, 10),
  show_rapid_points = TRUE
)

print(p_excfirst)
print(p_inhfirst)
