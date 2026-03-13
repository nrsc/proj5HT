suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(patchwork)
  library(scales)
})

## =========================================================
## Helpers (your existing ones; unchanged)
## =========================================================

add_smoothed_pct <- function(d,
                             k_smooth = 11,
                             smooth_method = c("rollmean","rollmedian"),
                             col_in = "pct_of_baseline") {
  smooth_method <- match.arg(smooth_method)
  d <- d %>% dplyr::arrange(time)
  y <- d[[col_in]]

  y_smooth <- switch(
    smooth_method,
    rollmean   = zoo::rollmean(y, k = k_smooth, fill = NA, align = "center"),
    rollmedian = zoo::rollmedian(y, k = k_smooth, fill = NA, align = "center")
  )

  d %>% dplyr::mutate(pct_smooth = as.numeric(y_smooth))
}

get_rapid_extrema <- function(d,
                              rapid_window = c(0, 10),
                              col = "pct_smooth") {
  d0 <- d %>% dplyr::filter(time >= rapid_window[1], time <= rapid_window[2])
  v <- d0[[col]]
  if (nrow(d0) == 0 || all(!is.finite(v))) {
    return(tibble::tibble(
      rapid_min_time = NA_real_, rapid_min_pct = NA_real_,
      rapid_max_time = NA_real_, rapid_max_pct = NA_real_
    ))
  }
  i_min <- which.min(v)
  i_max <- which.max(v)

  tibble::tibble(
    rapid_min_time = d0$time[i_min],
    rapid_min_pct  = d0[[col]][i_min],
    rapid_max_time = d0$time[i_max],
    rapid_max_pct  = d0[[col]][i_max]
  )
}

## =========================================================
## 1) Filter per-cell table to biphasic order of interest
## =========================================================

get_biphasic_cells <- function(per_cell2, biphasic_order_keep = c("inh_then_exc","exc_then_inh")) {
  stopifnot(all(c("cell_name","response_type","biphasic_order","assigned_subclass") %in% names(per_cell2)))

  per_cell2 %>%
    dplyr::filter(response_type == "biphasic") %>%
    dplyr::filter(biphasic_order %in% biphasic_order_keep) %>%
    dplyr::filter(!is.na(cell_name), !is.na(assigned_subclass))
}

## =========================================================
## 2) Build STACKED trace table (facet by subclass)
## =========================================================

prep_biphasic_stacked_df <- function(df_prepped,
                                     per_cell2,
                                     biphasic_order_keep = "inh_then_exc",
                                     xlim = c(-5, 20),
                                     k_smooth = 11,
                                     smooth_method = c("rollmean","rollmedian")) {
  smooth_method <- match.arg(smooth_method)

  stopifnot(all(c("cell_name","time","pct_of_baseline","assigned_subclass") %in% names(df_prepped)))

  keep_tbl <- get_biphasic_cells(per_cell2, biphasic_order_keep = biphasic_order_keep) %>%
    dplyr::select(cell_name, assigned_subclass) %>%
    dplyr::distinct()

  if (nrow(keep_tbl) == 0) stop("No biphasic cells found for biphasic_order_keep = ", biphasic_order_keep)

  d <- df_prepped %>%
    dplyr::semi_join(keep_tbl, by = "cell_name") %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(cell_name) %>%
    dplyr::arrange(time) %>%
    dplyr::group_modify(~add_smoothed_pct(.x, k_smooth = k_smooth, smooth_method = smooth_method)) %>%
    dplyr::ungroup()

  # Ensure subclass is present and consistent (df_prepped already has it; but keep_tbl is the truth)
  d <- d %>%
    dplyr::select(-assigned_subclass) %>%
    dplyr::left_join(keep_tbl, by = "cell_name")

  # trace index within each subclass (for stacking)
  ord <- d %>%
    dplyr::distinct(assigned_subclass, cell_name) %>%
    dplyr::arrange(assigned_subclass, cell_name) %>%
    dplyr::group_by(assigned_subclass) %>%
    dplyr::mutate(trace_idx = dplyr::row_number()) %>%
    dplyr::ungroup()

  d %>% dplyr::left_join(ord, by = c("assigned_subclass","cell_name"))
}

## =========================================================
## 3) Rapid markers (computed from SMOOTH; one min + one max)
## =========================================================

prep_rapid_markers_stacked <- function(d_stack,
                                       rapid_window = c(0, 10),
                                       stack_gap = 30) {
  # per cell, compute rapid extrema from pct_smooth
  marks <- d_stack %>%
    dplyr::group_by(assigned_subclass, cell_name, trace_idx) %>%
    dplyr::group_modify(~get_rapid_extrema(.x, rapid_window = rapid_window, col = "pct_smooth")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      offset = (trace_idx - 1) * stack_gap,
      y_min = rapid_min_pct + offset,
      y_max = rapid_max_pct + offset
    )

  marks
}

## =========================================================
## 4) STACKED plot by subclass (this replaces the bargraph)
## =========================================================

plot_biphasic_stacked_by_subclass <- function(df_prepped,
                                              per_cell2,
                                              biphasic_order_keep = "inh_then_exc",
                                              xlim = c(-5, 20),
                                              rapid_window = c(0, 10),
                                              k_smooth = 11,
                                              smooth_method = "rollmean",
                                              stack_gap = 30,
                                              show_raw = TRUE,
                                              show_smoothed = TRUE,
                                              show_rapid_points = TRUE) {

  d_stack <- prep_biphasic_stacked_df(
    df_prepped = df_prepped,
    per_cell2  = per_cell2,
    biphasic_order_keep = biphasic_order_keep,
    xlim = xlim,
    k_smooth = k_smooth,
    smooth_method = smooth_method
  ) %>%
    dplyr::mutate(
      offset = (trace_idx - 1) * stack_gap,
      y_raw    = pct_of_baseline + offset,
      y_smooth = pct_smooth + offset
    )

  # rapid markers (from smooth)
  marks <- NULL
  if (isTRUE(show_rapid_points)) {
    marks <- prep_rapid_markers_stacked(d_stack, rapid_window = rapid_window, stack_gap = stack_gap)
  }

  # color convention: RED = excitation (rapid max), BLUE = inhibition (rapid min)
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  base_lines <- d_stack %>%
    dplyr::distinct(assigned_subclass, cell_name, trace_idx, offset) %>%
    dplyr::mutate(y0 = 100 + offset)

  ggplot(d_stack, aes(x = time)) +
    facet_grid(assigned_subclass ~ ., scales = "free_y") +

    # rapid window shading
    annotate("rect",
             xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf,
             alpha = 0.06, fill = "grey30") +

    # per-trace baseline at 100 + offset
    geom_hline(data = base_lines, aes(yintercept = y0),
               linetype = 2, linewidth = 0.35, color = "grey45") +

    { if (isTRUE(show_raw))
      geom_line(aes(y = y_raw, group = cell_name),
                color = "grey75", linewidth = 0.25, alpha = 0.6)
      else NULL } +

    { if (isTRUE(show_smoothed))
      geom_line(aes(y = y_smooth, group = cell_name),
                color = "black", linewidth = 0.6, na.rm = TRUE)
      else NULL } +

    # rapid min/max points from smooth
    { if (!is.null(marks))
      geom_point(data = marks, aes(x = rapid_min_time, y = y_min),
                 color = col_inh, size = 1.7, na.rm = TRUE)
      else NULL } +
    { if (!is.null(marks))
      geom_point(data = marks, aes(x = rapid_max_time, y = y_max),
                 color = col_exc, size = 1.7, na.rm = TRUE)
      else NULL } +

    coord_cartesian(xlim = xlim) +
    labs(
      title = paste0("Rapid biphasic stacked traces (", biphasic_order_keep, ")"),
      subtitle = paste0("Window: [", xlim[1], ", ", xlim[2], "] s; rapid window: [",
                        rapid_window[1], ", ", rapid_window[2],
                        "] s; points = rapid min (blue) + rapid max (red)"),
      x = "Time (s)",
      y = "% baseline (100) + vertical offset"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = "grey70")
    )
}

## =========================================================
## 5) OPTIONAL: Zoomed patchwork (your old QC view) under the stacked plot
##     (kept, but NO bargraph anywhere)
## =========================================================

plot_zoom_cell_rapid_points <- function(df_prepped,
                                        cell_id,
                                        xlim = c(-5, 20),
                                        rapid_window = c(0, 10),
                                        k_smooth = 11,
                                        smooth_method = c("rollmean","rollmedian"),
                                        min_excursion_pct = 5,
                                        show_raw = TRUE,
                                        show_smoothed = TRUE) {
  smooth_method <- match.arg(smooth_method)

  d <- df_prepped %>%
    dplyr::filter(cell_name == cell_id) %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2])

  if (nrow(d) == 0) {
    return(ggplot() + theme_void() + labs(title = paste(cell_id, "(no data in xlim)")))
  }

  d <- add_smoothed_pct(d, k_smooth = k_smooth, smooth_method = smooth_method, col_in = "pct_of_baseline")
  rapid_pts <- get_rapid_extrema(d, rapid_window = rapid_window, col = "pct_smooth")

  col_exc <- "#d62728"
  col_inh <- "#1f77b4"
  thr_exc <- 100 + min_excursion_pct
  thr_inh <- 100 - min_excursion_pct

  ggplot(d, aes(x = time)) +
    annotate("rect", xmin = rapid_window[1], xmax = rapid_window[2],
             ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.4) +
    geom_hline(yintercept = thr_exc, linetype = 3, linewidth = 0.35, color = col_exc) +
    geom_hline(yintercept = thr_inh, linetype = 3, linewidth = 0.35, color = col_inh) +
    { if (show_raw)
      geom_line(aes(y = pct_of_baseline), linewidth = 0.25, alpha = 0.6, color = "grey55")
      else NULL } +
    { if (show_smoothed)
      geom_line(aes(y = pct_smooth), linewidth = 0.65, na.rm = TRUE, color = "black")
      else NULL } +
    { if (is.finite(rapid_pts$rapid_min_time) && is.finite(rapid_pts$rapid_min_pct))
      geom_point(data = rapid_pts, aes(x = rapid_min_time, y = rapid_min_pct),
                 color = col_inh, size = 1.9)
      else NULL } +
    { if (is.finite(rapid_pts$rapid_max_time) && is.finite(rapid_pts$rapid_max_pct))
      geom_point(data = rapid_pts, aes(x = rapid_max_time, y = rapid_max_pct),
                 color = col_exc, size = 1.9)
      else NULL } +
    coord_cartesian(xlim = xlim) +
    labs(title = cell_id, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 7)
    )
}

plot_patchwork_biphasic <- function(df_prepped,
                                    per_cell2,
                                    biphasic_order_keep = "inh_then_exc",
                                    xlim = c(-5, 20),
                                    rapid_window = c(0, 10),
                                    k_smooth = 11,
                                    smooth_method = c("rollmean","rollmedian"),
                                    min_excursion_pct = 5,
                                    n_show = 24,
                                    ncol = 6,
                                    seed = 1) {
  set.seed(seed)

  keep_tbl <- get_biphasic_cells(per_cell2, biphasic_order_keep = biphasic_order_keep)
  cell_ids <- keep_tbl$cell_name
  if (length(cell_ids) == 0) stop("No biphasic cells for order = ", biphasic_order_keep)
  if (length(cell_ids) > n_show) cell_ids <- sample(cell_ids, n_show)

  plots <- lapply(cell_ids, function(cid) {
    plot_zoom_cell_rapid_points(
      df_prepped, cid,
      xlim = xlim,
      rapid_window = rapid_window,
      k_smooth = k_smooth,
      smooth_method = smooth_method,
      min_excursion_pct = min_excursion_pct
    )
  })

  wrap_plots(plots, ncol = ncol) +
    plot_annotation(
      title = paste0("Rapid biphasic (", biphasic_order_keep, "): zoomed traces"),
      subtitle = paste0("Markers = rapid min (blue) + rapid max (red); xlim = [", xlim[1], ", ", xlim[2], "]")
    )
}

plot_biphasic_summary_no_bargraph <- function(df_prepped,
                                              per_cell2,
                                              biphasic_order_keep = "inh_then_exc",
                                              xlim = c(-5, 20),
                                              rapid_window = c(0, 10),
                                              k_smooth = 11,
                                              smooth_method = "rollmean",
                                              stack_gap = 30,
                                              n_show = 24,
                                              ncol_patch = 6) {

  p_stacked <- plot_biphasic_stacked_by_subclass(
    df_prepped = df_prepped,
    per_cell2  = per_cell2,
    biphasic_order_keep = biphasic_order_keep,
    xlim = xlim,
    rapid_window = rapid_window,
    k_smooth = k_smooth,
    smooth_method = smooth_method,
    stack_gap = stack_gap,
    show_raw = TRUE,
    show_smoothed = TRUE,
    show_rapid_points = TRUE
  )

  p_patch <- plot_patchwork_biphasic(
    df_prepped = df_prepped,
    per_cell2  = per_cell2,
    biphasic_order_keep = biphasic_order_keep,
    xlim = xlim,
    rapid_window = rapid_window,
    k_smooth = k_smooth,
    smooth_method = smooth_method,
    n_show = n_show,
    ncol = ncol_patch
  )

  p_stacked / p_patch + plot_layout(heights = c(2, 3))
}

## =========================================================
## Example usage
## =========================================================

# Stacked-only (what you asked for)
p_stacked_exc_inh <- plot_biphasic_stacked_by_subclass(
  df_prepped = df_prepped,
  per_cell2  = per_cell2,
  biphasic_order_keep = "exc_then_inh",
  xlim = c(-5, 30),
  rapid_window = c(0, 10),
  k_smooth = 3,
  smooth_method = "rollmean",
  stack_gap = 0
)
print(p_stacked_exc_inh)

p_stacked_inh_exc <- plot_biphasic_stacked_by_subclass(
  df_prepped = df_prepped,
  per_cell2  = per_cell2,
  biphasic_order_keep = "inh_then_exc",
  xlim = c(-5, 30),
  rapid_window = c(0, 10),
  k_smooth = 3,
  smooth_method = "rollmean",
  stack_gap = 0
)
print(p_stacked_inh_exc)


unique(per_cell2$rapid_pattern)
# #If you also want stacked + patchwork (still NO bargraph)
# p_combo_inh_exc <- plot_biphasic_summary_no_bargraph(
#   df_prepped = df_prepped,
#   per_cell2  = per_cell2,
#   biphasic_order_keep = "inh_then_exc",
#   xlim = c(-5, 30),
#   k_smooth = 3,
#   smooth_method = "rollmean",
#   stack_gap = 0,
#   n_show = 30,
#   ncol_patch = 6
# )
# print(p_combo_inh_exc)
# #If you also want stacked + patchwork (still NO bargraph)
# p_combo_exc_inh <- plot_biphasic_summary_no_bargraph(
#   df_prepped = df_prepped,
#   per_cell2  = per_cell2,
#   biphasic_order_keep = "exc_then_inh",
#   xlim = c(-5, 30),
#   k_smooth = 3,
#   smooth_method = "rollmean",
#   stack_gap = 0,
#   n_show = 30,
#   ncol_patch = 6
# )
# print(p_combo_exc_inh)
