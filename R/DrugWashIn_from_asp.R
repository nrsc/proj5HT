DrugWashIn_from_asp <- function(x,
                                mMD = NULL,
                                mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                                mMD_sheet = "Master_Macaque",
                                ttl_default = 5,
                                outlier_sd = 5,
                                bl_window = c(-5, 0),
                                plot_it = TRUE,
                                write_csv = TRUE,
                                csv_name = "-srtDrugWash.csv",
                                save_srt = FALSE,
                                verbose = TRUE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
  })

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

  # wash-in experiment labels (match your sheet usage)
  wash_exps <- c("Ket-wash-in", "WAYwash")

  # find testSpike protocol(s)
  meta_wash <- pick_asp_protocols(
    asp,
    rank = c("testSpike"),
    exp = wash_exps,
    stim_grepl = "spikeTTL|TTLstack"
  )

  if (nrow(meta_wash) == 0) {
    if (isTRUE(verbose)) message("DrugWashIn_from_asp: no testSpike wash-in protocol found.")
    return(NULL)
  }

  key_wash <- meta_wash$protocol_key[1]
  obj_wash <- asp$by_protocol[[key_wash]]

  # baseline protocol (same approach as spikePuff_from_asp)
  meta_bl <- pick_asp_protocols(
    asp,
    rank = c("1", "blSpike", "baseline", "bl", "Rank1"),
    stim_grepl = "spikeTTL|TTLstack"
  )
  if (nrow(meta_bl) == 0) meta_bl <- pick_asp_protocols(asp, stim_grepl = "spikeTTL")
  if (nrow(meta_bl) == 0) {
    if (isTRUE(verbose)) message("DrugWashIn_from_asp: no baseline protocol found; returning wash-only.")
  }

  key_bl <- if (nrow(meta_bl) > 0) meta_bl$protocol_key[1] else NA_character_
  obj_bl <- if (is.finite(match(key_bl, names(asp$by_protocol)))) asp$by_protocol[[key_bl]] else NULL

  wash_df <- obj_wash$spike_puff_output
  base_df <- if (!is.null(obj_bl)) obj_bl$spike_puff_output else NULL

  p <- NULL
  if (isTRUE(plot_it) && !is.null(base_df)) {
    df_plot <- dplyr::bind_rows(
      base_df %>% dplyr::select(time, percent_change) %>% dplyr::mutate(condition = "baseline (spikeTTL)"),
      wash_df %>% dplyr::select(time, percent_change) %>% dplyr::mutate(condition = "wash-in (testSpike)")
    )

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = percent_change, colour = condition)) +
      ggplot2::geom_hline(yintercept = 100, linetype = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::labs(
        title = paste0(asp$cell, " - ", as.character(obj_wash$pro_row$Exp)[1], " - ", as.character(obj_wash$pro_row$assigned_subclass)),
        x = "time (s)",
        y = "% of baseline"
      ) +
      ggplot2::theme_bw()

    print(p)
  }

  if (isTRUE(write_csv)) {
    rd <- if (!is.null(srt)) srt$rd else getwd()
    utils::write.csv(
      wash_df,
      file.path(rd, paste0(asp$cell, csv_name)),
      row.names = FALSE
    )
    if (isTRUE(verbose)) message("Saved: ", file.path(rd, paste0(asp$cell, csv_name)))
  }

  invisible(list(
    cell = asp$cell,
    baseline_key = key_bl,
    wash_key = key_wash,
    baseline = obj_bl,
    wash = obj_wash,
    plot = p
  ))
}
