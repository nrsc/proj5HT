run_asp_dfs_pipeline <- function(srt,
                                 mMD = NULL,
                                 mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                                 mMD_sheet = "Master_Macaque",
                                 ttl_default = 5,
                                 outlier_sd = 5,
                                 bl_window = c(-5, 0),
                                 save_srt = TRUE,
                                 run_spikePuff = TRUE,
                                 run_washin = TRUE,
                                 verbose = TRUE) {

  stop_if_not_srt(srt)

  # Ensure asp exists (heavy if missing, cheap if present)
  srt <- ensure_asp_in_srt(
    srt,
    mMD = mMD,
    mMD_path = mMD_path,
    mMD_sheet = mMD_sheet,
    ttl_default = ttl_default,
    outlier_sd = outlier_sd,
    bl_window = bl_window,
    save_srt = save_srt
  )

  if (is.null(srt$dfs$single) || !is.list(srt$dfs$single)) {
    srt$dfs$single <- list()
  }
  if (is.null(srt$dfs$drug_wash) || !is.list(srt$dfs$drug_wash)) {
    srt$dfs$drug_wash <- list()
  }

  # -----------------------------
  # spikePuff derived FROM asp
  # -----------------------------
  if (isTRUE(run_spikePuff)) {
    if (verbose) message("DFS pipeline: spikePuff_from_asp()")
    srt$dfs$single$spikePuff <- spikePuff_from_asp(
      srt,
      mMD = mMD,
      mMD_path = mMD_path,
      mMD_sheet = mMD_sheet,
      ttl_default = ttl_default,
      outlier_sd = outlier_sd,
      bl_window = bl_window,
      write_csv = TRUE,
      save_srt = save_srt,
      verbose = verbose
    )
  }

  # -----------------------------
  # Drug wash-in derived FROM asp
  # -----------------------------
  if (isTRUE(run_washin)) {
    if (verbose) message("DFS pipeline: DrugWashIn_from_asp()")
    srt$dfs$drug_wash$washin <- DrugWashIn_from_asp(
      srt,
      mMD = mMD,
      mMD_path = mMD_path,
      mMD_sheet = mMD_sheet,
      ttl_default = ttl_default,
      outlier_sd = outlier_sd,
      bl_window = bl_window,
      plot_it = FALSE,
      write_csv = TRUE,
      save_srt = save_srt,
      verbose = verbose
    )
  }

  if (isTRUE(save_srt)) {
    saveRDS(srt, file.path(srt$rd, paste0(srt$cell, "-srt.rds")))
  }

  srt$dfs
}
