run_asp_pipeline <- function(cell,
                             mMD = NULL,
                             mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                             mMD_sheet = "Master_Macaque",
                             ttl_default = 5,
                             outlier_sd = 5,
                             bl_window = c(-5, 0),
                             save_srt = TRUE,
                             make_concat_plot = TRUE,
                             make_single_plots = TRUE,
                             make_washin = TRUE,
                             verbose = TRUE) {

  got <- ensure_asp(
    cell,
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

  res <- list(asp = asp)

  if (isTRUE(make_concat_plot)) {
    res$concat <- plot_asp_concatenated(asp, verbose = verbose)
  }

  if (isTRUE(make_single_plots)) {
    keys <- names(asp$by_protocol)
    res$single <- lapply(keys, function(k) {
      plot_asp_single_protocol(asp, protocol_key = k, verbose = verbose)
    })
    names(res$single) <- keys
  }

  if (isTRUE(make_washin)) {
    res$washin <- DrugWashIn_from_asp(srt, mMD = mMD, plot_it = TRUE, write_csv = TRUE, verbose = verbose)
  }

  invisible(res)
}
