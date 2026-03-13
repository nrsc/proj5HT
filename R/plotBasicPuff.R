plotBasicPuff <- function(x,
                          show_plot = TRUE,
                          save_plot = TRUE,
                          bin_seconds = 1,
                          xlim = c(-20, 50),
                          figs_root = "figs/basicPuff",
                          colours = c(
                            percent_change = "red",
                            instRate = "blue"
                          )) {

  # ---- helpers ----
  `%||%` <- function(a, b) if (!is.null(a)) a else b

  slugify <- function(s) {
    s <- as.character(s)
    s <- s[!is.na(s) & nzchar(s)]
    if (!length(s)) return("unknown")
    s <- s[1]
    s <- gsub("\\s+", "_", s)
    s <- gsub("[\\[\\]\\{\\}()]", "", s)
    s <- gsub("[/\\\\:;,*?\"<>|]+", "-", s)
    s <- gsub("[^A-Za-z0-9._+-]+", "_", s)
    s <- gsub("_+", "_", s)
    s <- gsub("(^_+|_+$)", "", s)
    if (!nzchar(s)) "unknown" else s
  }

  # ---- checks ----
  if (is.null(x$dfs$spikeTTL) || is.null(x$dfs$spikeTTL$spike_puff_output)) {
    message("No spikeTTL to plot.")
    return(NULL)
  }

  df <- x$dfs$spikeTTL$spike_puff_output

  # ---- title parts ----
  cell <- x$md$cell_name %||% x$cell %||% "unknown_cell"

  assigned <- unique(na.omit(df$assigned_subclass))
  assigned <- if (length(assigned)) paste(assigned, collapse = ",") else "unknown_subclass"

  puff <- unique(na.omit(df$puff))
  puff <- if (length(puff)) paste(puff, collapse = ",") else "unknown_puff"

  bath <- if ("bath" %in% names(df)) unique(na.omit(df$bath)) else NULL
  bath <- if (!is.null(bath) && length(bath)) paste(bath, collapse = ",") else NULL

  title_string <- if (!is.null(bath)) {
    paste(cell, assigned, puff, paste0("Bath:", bath), sep = " - ")
  } else {
    paste(cell, assigned, puff, sep = " - ")
  }

  # ---- bin to seconds ----
  df_binned <- df |>
    dplyr::mutate(time_bin = floor(time / bin_seconds) * bin_seconds) |>
    dplyr::group_by(time_bin) |>
    dplyr::summarise(
      percent_change = mean(percent_change, na.rm = TRUE),
      instRate = mean(instRate, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- scale factor for dual axis ----
  pc_range <- range(df_binned$percent_change, na.rm = TRUE)
  ir_range <- range(df_binned$instRate, na.rm = TRUE)

  scale_factor <- diff(pc_range) / diff(ir_range)

  # ---- plot ----
  gg <- ggplot2::ggplot(df_binned, ggplot2::aes(x = time_bin)) +
    # % change (left axis)
    ggplot2::geom_line(
      ggplot2::aes(y = percent_change),
      colour = colours["percent_change"],
      linewidth = 1
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = percent_change),
      colour = colours["percent_change"],
      size = 1.5
    ) +
    # instRate (right axis, scaled)
    ggplot2::geom_line(
      ggplot2::aes(y = instRate * scale_factor),
      colour = colours["instRate"],
      linewidth = 1
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_hline(yintercept = 100, linetype = 3) +
    ggplot2::scale_y_continuous(
      name = "% Change (vs baseline mean)",
      sec.axis = ggplot2::sec_axis(
        ~ . / scale_factor,
        name = "Instantaneous firing rate (Hz)"
      )
    ) +
    ggplot2::coord_cartesian(xlim = xlim) +
    ggplot2::labs(
      title = title_string,
      x = "Time (s)"
    ) +
    ggplot2::theme_bw()

  if (isTRUE(show_plot)) print(gg)

  if (isTRUE(save_plot)) {
    fname <- paste0(slugify(title_string), ".png")

    # project figs: per-puff folder
    puff_dir <- slugify(puff)
    rel_dir <- file.path(figs_root, puff_dir)
    if (!dir.exists(rel_dir)) dir.create(rel_dir, recursive = TRUE)

    ggplot2::ggsave(file.path(rel_dir, fname), plot = gg, width = 6, height = 5)

    # x$rd: flat only
    if (!is.null(x$rd) && dir.exists(x$rd)) {
      ggplot2::ggsave(file.path(x$rd, fname), plot = gg, width = 6, height = 5)
    }
  }
  invisible(gg)
}
