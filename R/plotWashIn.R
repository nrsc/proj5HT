plotWashIn <- function(x,
                       show_plot = TRUE,
                       save_plot = TRUE,
                       bin_seconds = 1,
                       xlim = c(-20, 60),
                       figs_root = "figs/washIn",
                       colours = c(
                         "baseline" = "grey40",
                         "wash-in" = "red"
                       )) {

  # --- helpers ---
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  slugify <- function(s) {
    s <- as.character(s)
    s <- s[!is.na(s) & nzchar(s)]
    if (!length(s)) return("unknown")
    s <- s[1]
    s <- gsub("\\s+", "_", s)                  # spaces -> _
    s <- gsub("[\\[\\]\\{\\}()]", "", s)       # remove brackets safely
    s <- gsub("[/\\\\:;,*?\"<>|]+", "-", s)    # forbidden chars -> -
    s <- gsub("[^A-Za-z0-9._+-]+", "_", s)     # everything else -> _
    s <- gsub("_+", "_", s)                    # collapse ___
    s <- gsub("(^_+|_+$)", "", s)              # trim _
    if (!nzchar(s)) "unknown" else s
  }


  # --- checks ---
  if (is.null(x$dfs$DrugWash) || is.null(x$dfs$DrugWash$spike_puff_output)) {
    message("No drug wash-in to plot.")
    return(NULL)
  }
  if (is.null(x$dfs$spikeTTL) || is.null(x$dfs$spikeTTL$spike_puff_output)) {
    message("No baseline spikeTTL to plot.")
    return(NULL)
  }

  baseline_df <- x$dfs$spikeTTL$spike_puff_output
  wash_df     <- x$dfs$DrugWash$spike_puff_output

  # --- title parts ---
  cell <- x$md$cell_name %||% x$cell %||% "unknown_cell"

  assigned <- unique(na.omit(c(baseline_df$assigned_subclass, wash_df$assigned_subclass)))
  assigned <- if (length(assigned)) paste(assigned, collapse = ",") else "unknown_subclass"

  puff <- unique(na.omit(c(baseline_df$puff, wash_df$puff)))
  puff <- if (length(puff)) paste(puff, collapse = ",") else "unknown_puff"

  bath <- unique(na.omit(wash_df$bath))
  bath <- if (length(bath)) paste(bath, collapse = ",") else "unknown_bath"

  exp_name <- unique(na.omit(wash_df$expCon))
  exp_name <- if (length(exp_name)) exp_name[1] else "WashIn"

  title_string <- paste(cell, assigned, puff, paste0("WashIn:", bath), sep = " - ")

  # --- build tidy df for plotting ---
  df_plot <- dplyr::bind_rows(
    baseline_df |>
      dplyr::transmute(time, instRate, percent_change,
                       condition = "baseline"),
    wash_df |>
      dplyr::transmute(time, instRate, percent_change,
                       condition = "wash-in")
  )

  # bin to seconds
  df_plot <- df_plot |>
    dplyr::mutate(time_bin = floor(time / bin_seconds) * bin_seconds) |>
    dplyr::group_by(time_bin, condition) |>
    dplyr::summarise(
      instRate = mean(instRate, na.rm = TRUE),
      percent_change = mean(percent_change, na.rm = TRUE),
      .groups = "drop"
    )

  # long format for two panels
  df_long <- tidyr::pivot_longer(
    df_plot,
    cols = c(instRate, percent_change),
    names_to = "metric",
    values_to = "value"
  )

  df_long$metric <- factor(
    df_long$metric,
    levels = c("instRate", "percent_change"),
    labels = c("Instantaneous firing rate (Hz)", "% Change (vs baseline mean)")
  )

  # plot
  gg <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(x = time_bin, y = value, colour = condition, group = condition)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(
      data = subset(df_long, metric == "% Change (vs baseline mean)"),
      ggplot2::aes(yintercept = 100),
      linetype = 3,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
    ggplot2::geom_point(size = 1.5, na.rm = TRUE) +
    ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(title = title_string, x = "Time (s)", y = NULL, colour = NULL) +
    ggplot2::theme_bw()

  if (isTRUE(show_plot)) print(gg)

  if (isTRUE(save_plot)) {
    # folder structure: <figs_root>/<puff>_<bath>/
    condition_dir <- slugify(paste(puff, bath, sep = "_"))

    rel_dir <- file.path(figs_root, condition_dir)
    if (!dir.exists(rel_dir)) dir.create(rel_dir, recursive = TRUE)

    fname <- paste0(slugify(title_string), ".png")

    ggplot2::ggsave(file.path(rel_dir, fname), plot = gg, width = 6, height = 8)

    if (!is.null(x$rd) && dir.exists(x$rd)) {
      ggplot2::ggsave(file.path(x$rd, fname), plot = gg, width = 6, height = 8)
    }

  }

  invisible(gg)
}
