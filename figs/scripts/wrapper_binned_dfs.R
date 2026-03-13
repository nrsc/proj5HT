plot_binned_cell_traces <- function(df,
                                    x_col = "time",
                                    y_col = "percent_change",
                                    cell_col = "cell_name",
                                    color_col = "assigned_subclass",
                                    facet_col = "assigned_subclass",
                                    extra_group_cols = c("Species"),
                                    bin_width = 2,
                                    use_existing_x_bin = FALSE,
                                    x_bin_col = "x_bin",
                                    summary_fun = mean,
                                    na_rm = TRUE,
                                    xlim = c(-10, 50),
                                    ylim = c(0, 300),
                                    title = "5HT response across subclass",
                                    x_lab = "Time (s)",
                                    y_lab = "Percent Change",
                                    line_alpha = 0.8,
                                    line_size = 0.5,
                                    show_legend = FALSE) {

  # basic column check
  cols_needed <- c(x_col, y_col, cell_col, color_col, facet_col, extra_group_cols)
  cols_needed <- unique(cols_needed[!is.na(cols_needed) & nzchar(cols_needed)])
  missing_cols <- setdiff(cols_needed, colnames(df))

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # make / use x_bin
  d0 <- df

  if (use_existing_x_bin) {
    if (!x_bin_col %in% colnames(d0)) {
      stop("use_existing_x_bin = TRUE, but `", x_bin_col, "` is not in df")
    }
    d0[[".plot_x_bin"]] <- d0[[x_bin_col]]
  } else {
    d0[[".plot_x_bin"]] <- floor(d0[[x_col]] / bin_width) * bin_width
  }

  # grouping columns for summary
  group_cols <- c(".plot_x_bin", cell_col, color_col, extra_group_cols)
  group_cols <- unique(group_cols[!is.na(group_cols) & nzchar(group_cols)])

  # ensure facet column survives summarise if needed
  if (!facet_col %in% group_cols) {
    group_cols <- c(group_cols, facet_col)
  }

  # summarise
  d_sum <- d0 %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      y = summary_fun(.data[[y_col]], na.rm = na_rm),
      .groups = "drop"
    )

  # plot
  p <- ggplot2::ggplot(
    d_sum,
    ggplot2::aes(
      x = .data[[".plot_x_bin"]],
      y = .data[["y"]],
      group = .data[[cell_col]],
      colour = .data[[color_col]]
    )
  ) +
    ggplot2::geom_line(alpha = line_alpha, linewidth = line_size) +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab,
      colour = color_col
    ) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::theme_minimal()

  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)))
  }

  if (!show_legend) {
    p <- p + ggplot2::theme(legend.position = "none")
  }

  return(p)
}
