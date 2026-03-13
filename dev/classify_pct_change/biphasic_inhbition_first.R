suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(zoo)
  library(patchwork)
  library(scales)
})

## =========================================================
## Helpers
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
## 1) Filter per-cell table to inhibition->excitation biphasic
## =========================================================

get_inh_then_exc_cells <- function(per_cell2) {
  stopifnot(all(c("cell_name","response_type","biphasic_order","assigned_subclass") %in% names(per_cell2)))

  per_cell2 %>%
    dplyr::filter(response_type == "biphasic") %>%
    dplyr::filter(biphasic_order == "inh_then_exc") %>%
    dplyr::filter(!is.na(cell_name))
}

## =========================================================
## 2) Zoomed trace plot (single cell), only rapid points
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
  if (!"pct_of_baseline" %in% names(d)) {
    stop("df_prepped must include pct_of_baseline (run your prep that makes 100*instRate/baseline).")
  }

  d <- add_smoothed_pct(d, k_smooth = k_smooth, smooth_method = smooth_method, col_in = "pct_of_baseline")
  rapid_pts <- get_rapid_extrema(d, rapid_window = rapid_window, col = "pct_smooth")

  # Your convention: RED = excitation, BLUE = inhibition
  col_exc <- "#d62728"
  col_inh <- "#1f77b4"

  # Thresholds for “real excursion”
  thr_exc <- 100 + min_excursion_pct
  thr_inh <- 100 - min_excursion_pct

  # Title line (short for patchwork)
  title_txt <- cell_id

  p <- ggplot(d, aes(x = time)) +
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
    # Only rapid biphasic markers:
    { if (is.finite(rapid_pts$rapid_min_time) && is.finite(rapid_pts$rapid_min_pct))
      geom_point(data = rapid_pts, aes(x = rapid_min_time, y = rapid_min_pct),
                 color = col_inh, size = 1.9)
      else NULL } +
    { if (is.finite(rapid_pts$rapid_max_time) && is.finite(rapid_pts$rapid_max_pct))
      geom_point(data = rapid_pts, aes(x = rapid_max_time, y = rapid_max_pct),
                 color = col_exc, size = 1.9)
      else NULL } +
    coord_cartesian(xlim = xlim) +
    labs(title = title_txt, x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 7)
    )

  p
}

## =========================================================
## 3) Patchwork of many inhibition-first biphasic traces (zoomed)
## =========================================================

plot_patchwork_inh_then_exc <- function(df_prepped,
                                        per_cell2,
                                        xlim = c(-5, 20),
                                        rapid_window = c(0, 10),
                                        k_smooth = 11,
                                        smooth_method = c("rollmean","rollmedian"),
                                        min_excursion_pct = 5,
                                        n_show = 24,
                                        ncol = 6,
                                        seed = 1) {
  set.seed(seed)

  inh_first <- get_inh_then_exc_cells(per_cell2)
  if (nrow(inh_first) == 0) stop("No inhibition->excitation biphasic cells found (check biphasic_order labels).")

  # choose cells (random sample for QC browsing)
  cell_ids <- inh_first$cell_name
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
      title = "Rapid biphasic (inh → exc): zoomed traces",
      subtitle = paste0("Window: [", xlim[1], ", ", xlim[2], "] s; rapid window: [",
                        rapid_window[1], ", ", rapid_window[2],
                        "] s; markers = rapid min (blue) + rapid max (red)")
    )
}

## =========================================================
## 4) Stacked composition figure for inhibition-first biphasic
##     (how many of these cells per subclass)
## =========================================================

plot_stacked_inh_first_by_subclass <- function(per_cell2,
                                               min_n_label = 3) {
  inh_first <- get_inh_then_exc_cells(per_cell2)

  if (nrow(inh_first) == 0) stop("No inhibition->excitation biphasic cells found.")

  # Composition of inhibition-first biphasic cells across subclasses
  # Plot as proportions + show n in labels (useful for talks/QC)
  counts <- inh_first %>%
    dplyr::count(assigned_subclass, name = "n") %>%
    dplyr::mutate(frac = n / sum(n))

  ggplot(counts, aes(x = assigned_subclass, y = frac)) +
    geom_col(width = 0.8) +
    geom_text(aes(label = ifelse(n >= min_n_label, paste0("n=", n), "")),
              vjust = -0.4, size = 3) +
    scale_y_continuous(labels = percent_format(), limits = c(0, max(counts$frac) * 1.15)) +
    labs(
      x = NULL,
      y = "Proportion of inh→exc biphasic cells",
      title = "Inhibition-first biphasic cells: subclass composition"
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

## =========================================================
## 5) Combined “stacked + patchwork” layout
## =========================================================

plot_inh_first_summary <- function(df_prepped,
                                   per_cell2,
                                   xlim = c(-5, 20),
                                   n_show = 24,
                                   ncol_patch = 6,
                                   k_smooth = 11,
                                   smooth_method = "rollmean",
                                   min_excursion_pct = 5) {
  p_stack <- plot_stacked_inh_first_by_subclass(per_cell2)
  p_patch <- plot_patchwork_inh_then_exc(df_prepped, per_cell2,
                                         xlim = xlim,
                                         n_show = n_show,
                                         ncol = ncol_patch,
                                         k_smooth = k_smooth,
                                         smooth_method = smooth_method,
                                         min_excursion_pct = min_excursion_pct)

  p_stack / p_patch + plot_layout(heights = c(1, 3))
}

## =========================================================
## Example usage
## =========================================================

# 1) Make the combined figure
p_inh_first <- plot_inh_first_summary(
  df_prepped = df_prepped,
  per_cell2  = per_cell2,
  xlim = c(-5, 20),
  n_show = 30,          # how many traces to show
  ncol_patch = 6,
  k_smooth = 3,
  smooth_method = "rollmean",
  min_excursion_pct = 5
)

print(p_inh_first)

# 2) If you want *all* inh→exc traces saved as pages (QC):
#    (Optional: requires ggplot2 + patchwork; uses ggsave)
# save_inh_first_pages <- function(df_prepped, per_cell2,
#                                  out_dir = "qc_inh_then_exc_pages",
#                                  page_size = 24, ncol = 6,
#                                  xlim = c(-5, 20),
#                                  ...) {
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#   inh_first <- get_inh_then_exc_cells(per_cell2)
#   ids <- inh_first$cell_name
#   n_pages <- ceiling(length(ids) / page_size)
#
#   for (pg in seq_len(n_pages)) {
#     idx <- ((pg - 1) * page_size + 1):min(pg * page_size, length(ids))
#     p <- plot_patchwork_inh_then_exc(df_prepped, per_cell2 %>% dplyr::filter(cell_name %in% ids[idx]),
#                                      xlim = xlim, n_show = length(idx), ncol = ncol, ...)
#     fn <- file.path(out_dir, sprintf("inh_then_exc_page_%02d.png", pg))
#     ggsave(fn, p, width = 14, height = 8, dpi = 200)
#   }
#   message("Saved ", n_pages, " pages to: ", out_dir)
# }

#save_inh_first_pages(df_prepped, per_cell2, k_smooth = 2, smooth_method = "rollmean", min_excursion_pct = 15)
