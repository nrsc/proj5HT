spikePuff_from_asp <- function(x,
                               mMD = NULL,
                               mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                               mMD_sheet = "Master_Macaque",
                               ttl_default = 5,
                               outlier_sd = 5,
                               bl_window = c(-5, 0),
                               write_csv = TRUE,
                               csv_name = "-srtPuff.csv",
                               save_srt = FALSE,
                               verbose = TRUE,
                               # NEW
                               plot_it = FALSE,
                               plot_mode = c("percent_change", "log2_fc"),
                               min_rate_hz = 0.05) {

  suppressPackageStartupMessages({
    library(dplyr)
  })

  plot_mode <- match.arg(plot_mode)

  got <- ensure_asp(
    x,
    mMD = mMD,
    mMD_path = mMD_path,
    mMD_sheet = mMD_sheet,
    ttl_default = ttl_default,
    outlier_sd = outlier_sd,
    bl_window = bl_window,
    save_srt = save_srt
  )
  srt <- got$srt
  asp <- got$asp

  # choose a baseline-ish spikeTTL protocol row (mirror your old logic)
  meta <- pick_asp_protocols(
    asp,
    rank = c("1", "blSpike", "baseline", "bl", "Rank1"),
    stim_grepl = "spikeTTL|TTLstack"
  )

  if (nrow(meta) == 0) {
    meta <- pick_asp_protocols(asp, stim_grepl = "spikeTTL")
  }
  if (nrow(meta) == 0) {
    if (isTRUE(verbose)) message("spikePuff_from_asp: no suitable spikeTTL protocol found.")
    return(NULL)
  }

  key <- meta$protocol_key[1]
  obj <- asp$by_protocol[[key]]

  # main output df
  df <- obj$spike_puff_output

  # optional simple plot
  p <- NULL
  if (isTRUE(plot_it)) {
    suppressPackageStartupMessages({
      library(ggplot2)
    })

    # pick y
    if (plot_mode == "percent_change") {
      y <- df$percent_change
      ylab <- "% of baseline"
      hline <- 100
    } else {
      # log2_fc from instRate/blMean (stable with min_rate_hz)
      bl <- obj$blMean
      y <- log2(pmax(df$instRate, min_rate_hz) / bl)
      ylab <- "log2(instRate / baseline)"
      hline <- 0
    }

    dfp <- df %>%
      dplyr::mutate(.y = y) %>%
      dplyr::filter(is.finite(.data$time), is.finite(.data$.y)) %>%
      dplyr::arrange(.data$time)

    p <- ggplot2::ggplot(dfp, ggplot2::aes(x = time, y = .y)) +
      ggplot2::geom_hline(yintercept = hline, linewidth = 0.3, linetype = 3) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, linetype = 2) +
      ggplot2::geom_line(linewidth = 0.4, na.rm = TRUE) +
      ggplot2::labs(
        title = paste0(asp$cell, " | spikePuff_from_asp"),
        subtitle = key,
        x = "time (s)",
        y = ylab
      ) +
      ggplot2::theme_classic()

    print(p)
  }

  out <- list(
    cell = asp$cell,
    protocol_key = key,
    pro_row = obj$pro_row,
    spike_puff_output = df,
    raw_trace = obj$raw_trace,
    blMean = obj$blMean,
    plot = p
  )

  if (isTRUE(write_csv)) {
    rd <- if (!is.null(srt)) srt$rd else getwd()
    utils::write.csv(
      out$spike_puff_output,
      file.path(rd, paste0(asp$cell, csv_name)),
      row.names = FALSE
    )
    if (isTRUE(verbose)) message("Saved: ", file.path(rd, paste0(asp$cell, csv_name)))
  }

  out
}
