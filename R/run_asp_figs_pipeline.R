run_asp_figs_pipeline <- function(srt,
                                  make_concat_plot = TRUE,
                                  make_single_plots = TRUE,
                                  concat_args = list(),
                                  single_args = list(),
                                  verbose = TRUE) {
  stop_if_not_srt(srt)

  if (is.null(srt$figs) || !is.list(srt$figs)) srt$figs <- list()
  if (is.null(srt$figs$asp) || !is.list(srt$figs$asp)) srt$figs$asp <- list()

  asp <- srt$dfs$asp
  if (is.null(asp) || is.null(asp$by_protocol) || length(asp$by_protocol) == 0) {
    stop("FIGS pipeline: srt$dfs$asp is missing/empty. Run run_asp_dfs_pipeline() first.")
  }

  # Concatenated plot
  if (isTRUE(make_concat_plot)) {
    if (verbose) message("FIGS pipeline: plot_asp_concatenated()")
    srt$figs$asp$concat <- do.call(
      plot_asp_concatenated,
      c(list(asp = asp, verbose = verbose), concat_args)
    )
  }

  # Single protocol plots
  if (isTRUE(make_single_plots)) {
    keys <- names(asp$by_protocol)
    if (verbose) message("FIGS pipeline: plot_asp_single_protocol() x ", length(keys))

    out <- lapply(seq_along(keys), function(i) {
      k <- keys[i]
      do.call(
        plot_asp_single_protocol,
        c(list(asp = asp, protocol_key = k, verbose = verbose), single_args)
      )
    })
    names(out) <- keys
    srt$figs$asp$single <- out
  }

  srt$figs
}
