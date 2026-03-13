suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

plot_cell_trace <- function(df, per_cell, cell_id,
                            baseline_value = 100,
                            protocol_col = "protocol",
                            time_col = "time",
                            y_col = "percent_change",
                            t0 = 0) {

  stopifnot(all(c(cell_id) %in% df$cell_name))
  stopifnot(all(c(protocol_col, time_col, y_col) %in% names(df)))
  stopifnot(all(c("cell_name", "response_call") %in% names(per_cell)))

  dfi <- df %>%
    dplyr::filter(cell_name == cell_id) %>%
    dplyr::mutate(
      protocol_std = tolower(trimws(.data[[protocol_col]])),
      is_baseline = protocol_std == "baseline"
    )

  call_i <- per_cell %>%
    dplyr::filter(cell_name == cell_id) %>%
    dplyr::distinct(response_call) %>%
    dplyr::pull(response_call)

  call_i <- ifelse(length(call_i) == 0, "NA", call_i[1])

  ggplot(dfi, aes(x = .data[[time_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = baseline_value, linetype = "dashed") +
    geom_vline(xintercept = t0, linetype = "dotted") +
    geom_point(aes(shape = is_baseline), alpha = 0.6, size = 1.4) +
    geom_line(alpha = 0.4) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
    xlim(-5,NA) +
    labs(
      title = paste0(cell_id, "  |  call = ", call_i),
      x = "Time (s)",
      y = "percent_change (baseline = 100)"
    ) +
    theme_bw()
}

# ---- Loop: make a bunch of plots ----
# Choose a set of cells to review:
cells_to_check <- res$per_cell %>%
  dplyr::filter(!is.na(response_call)) %>%
  dplyr::arrange(response_call) %>%
  dplyr::pull(cell_name) %>%
  unique()

# Optionally only check N cells:
cells_to_check <- head(cells_to_check, 10)

# Print sequentially
for (cid in cells_to_check) {
  print(plot_cell_trace(df, res$per_cell, cid))
}

# ---- Save to a folder (recommended) ----
out_dir <- "qc_traces"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (cid in cells_to_check) {
  p <- plot_cell_trace(df, res$per_cell, cid)
  ggsave(
    filename = file.path(out_dir, paste0(cid, "_qc.png")),
    plot = p,
    width = 9, height = 4, dpi = 150
  )
}
