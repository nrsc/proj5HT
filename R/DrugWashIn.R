#' analysis function for the ketanserin / WAY wash-in
#'
#' @param x srt object (from loadHCT/nnest)
#' @param mMD optional preloaded Master_UI_data sheet (recommended for loops)
#' @param mMD_path path to ODS if mMD not provided
#' @param mMD_sheet sheet name
#' @param ttl_default default TTL time if missing
#' @param outlier_sd remove points > outlier_sd SD from mean instRate
#' @param bl_window baseline window (seconds) used for blMean, e.g. c(-5, 0)
#' @param plot_it logical; make plots
#' @param write_csv logical; write csv to srt$rd
#' @return list (or NULL if not a washin cell / no usable data)
#' @export
DrugWashIn <- function(x,
                       mMD = NULL,
                       mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                       mMD_sheet = "Master_Macaque",
                       ttl_default = 5,
                       outlier_sd = 5,
                       bl_window = c(-5, 0),
                       plot_it = TRUE,
                       write_csv = TRUE) {

  stopifnot(is.list(x))

  MD <- projHCT::sheets$MD

  # Preload master sheet once if not provided
  if (is.null(mMD)) {
    mMD <- readODS::read_ods(mMD_path, sheet = mMD_sheet)
    mMD <- mMD[!is.na(mMD$cell_name), , drop = FALSE]
  }

  unique(mMD$cell_name)

  # Determine cell name robustly
  cell <- NULL
  if (!is.null(x$cell)) cell <- x$cell
  if (is.null(cell) && !is.null(x$md$cell)) cell <- x$md$cell
  if (is.null(cell) && !is.null(x$md$cell_name)) cell <- x$md$cell_name
  if (is.null(cell)) stop("Cannot determine cell name from x (srt).")

  # Wash-in cell list
  KWI_cells <- unique(mMD$cell_name[mMD$Exp == "Ket-wash-in"])
  WWI_cells <- unique(mMD$cell_name[mMD$Exp == "WAYwash"])
  washin_cells <- unique(c(KWI_cells, WWI_cells))

  if (!(cell %in% washin_cells)) {
    message("DrugWashIn: cell ", cell, " is not a Ket-wash-in or WAYwash experiment — skipping.")
    return(NULL)
  }

  # Require spikeTTL (baseline) if you want to compare/plot alongside
  spikeTTL <- x$dfs$spikeTTL
  # (not strictly required for analysis, but your original plots use it)
  if (is.null(spikeTTL)) message("Note: x$dfs$spikeTTL is NULL (baseline plot will be skipped).")

  # Match MD row
  md_idx <- which(grepl(cell, MD$cell_name))
  if (length(md_idx) == 0) {
    message("No MD match for cell: ", cell)
    return(NULL)
  }
  if (length(md_idx) > 1) {
    message("Multiple MD matches for cell ", cell, " (using first).")
    md_idx <- md_idx[1]
  }
  iMD <- MD[md_idx, , drop = FALSE]

  # Match Master row(s) for this cell
  m_idx <- which(grepl(cell, mMD$cell_name))
  if (length(m_idx) == 0) {
    message("No Master_UI_data match for cell: ", cell)
    return(NULL)
  }
  iMerged <- mMD[m_idx, , drop = FALSE]

  # Pick testSpike row for wash-in analysis
  if (!("Rank" %in% names(iMerged))) {
    message("Missing 'Rank' column in Master sheet for cell: ", cell)
    return(NULL)
  }
  test_idx <- which(iMerged$Rank == "testSpike")
  if (length(test_idx) == 0) {
    message("No Rank == 'testSpike' row for wash-in cell: ", cell)
    return(NULL)
  }
  iM0 <- iMerged[test_idx[1], , drop = FALSE]

  # TTL time
  ttlTime <- iM0$TTLtime
  if (is.na(ttlTime)) ttlTime <- ttl_default

  # Parse sweeps and protocol
  pro <- as.character(iM0$stimulus_description)
  spike_puff_sweeps <- as.character(iM0$wrapped_sweeps)
  spike_puff_sweeps <- gsub("\\[|\\]", "", spike_puff_sweeps)

  sp <- strsplit(spike_puff_sweeps, ",")[[1]]
  sp <- trimws(sp)
  sN <- as.numeric(sp) + 1

  # Extract NWB sweeps
  nwb <- nphys::extractNWB(iMD$nwb_file, sweeps = sN)
  if (is.null(names(nwb$sweeps))) nwb$sweeps <- as.data.frame(nwb$sweeps)

  sweep_names <- names(nwb$sweeps)

  # Durations (seconds)
  dur <- vapply(sweep_names, function(nm) {
    sw <- nwb$sweeps[[nm]]
    sw <- sw[!is.na(sw)]
    sr <- nwb$sampling_rate[[nm]]
    tail(nphys::convertIndex(sw, sr), 1)
  }, numeric(1))

  # APanalysis per sweep
  sweep_peaks <- lapply(sweep_names, function(nm) {
    sweep <- nwb$sweeps[[nm]]
    sweep <- sweep[!is.na(sweep)]
    sr <- nwb$sampling_rate[[nm]]
    nphys::APanalysis(sweep, rate = sr)
  })
  names(sweep_peaks) <- sweep_names

  # Name sweeps by protocol chunks
  pro_parts <- unlist(strsplit(pro, "\\+"))
  if (length(pro_parts) != length(sweep_peaks)) {
    names(sweep_peaks) <- c(pro_parts, rep("trail", length(sweep_peaks) - length(pro_parts)))
  } else {
    names(sweep_peaks) <- pro_parts
  }
  names(dur) <- names(sweep_peaks)

  # Build aligned dataframes
  dfBl <- dfTTL <- dfTrail <- NULL

  if ("baseline" %in% names(sweep_peaks)) {
    fp <- as.data.frame(sweep_peaks[[which(names(sweep_peaks) == "baseline")[1]]])
    fp <- fp[-nrow(fp), , drop = FALSE]
    base_dur <- as.numeric(dur[which(names(dur) == "baseline")[1]])
    dfBl <- data.frame(
      time = fp$epoch_t - (base_dur + ttlTime),
      instRate = fp$instRate
    )
  }

  if ("spikeTTL" %in% names(sweep_peaks)) {
    fp <- as.data.frame(sweep_peaks[[which(names(sweep_peaks) == "spikeTTL")[1]]])
    fp <- fp[-nrow(fp), , drop = FALSE]
    dfTTL <- data.frame(
      time = fp$epoch_t - ttlTime,
      instRate = fp$instRate
    )
  } else {
    message("No spikeTTL segment found in this protocol for cell: ", cell)
    return(NULL)
  }

  if ("trail" %in% names(sweep_peaks)) {
    trail_idx <- which(names(sweep_peaks) == "trail")

    if (length(trail_idx) >= 2) {
      fps <- lapply(trail_idx, function(ii) {
        d <- as.data.frame(sweep_peaks[[ii]])
        d <- d[-nrow(d), , drop = FALSE]
        d
      })

      pre_idx <- pmax(trail_idx - 1, 1)
      pre_dur <- as.numeric(dur[pre_idx])

      lens <- vapply(fps, nrow, integer(1))
      which_trail <- rep(seq_along(fps), times = lens)

      fp_all <- dplyr::bind_rows(fps)
      shift <- pre_dur[which_trail] * which_trail
      fp_all$epoch_t <- fp_all$epoch_t + shift - ttlTime

      dfTrail <- data.frame(time = fp_all$epoch_t, instRate = fp_all$instRate)
    } else {
      fp <- as.data.frame(sweep_peaks[[trail_idx[1]]])
      fp <- fp[-nrow(fp), , drop = FALSE]
      sp_dur <- as.numeric(dur[which(names(dur) == "spikeTTL")[1]])
      fp$epoch_t <- fp$epoch_t + sp_dur - ttlTime
      dfTrail <- data.frame(time = fp$epoch_t, instRate = fp$instRate)
    }
  }

  # Combine
  tst <- dfTTL
  if (!is.null(dfBl))    tst <- rbind(dfBl, tst)
  if (!is.null(dfTrail)) tst <- rbind(tst, dfTrail)

  # Outlier removal
  mu <- mean(tst$instRate, na.rm = TRUE)
  sig <- sd(tst$instRate, na.rm = TRUE)
  if (is.finite(sig) && sig > 0) {
    keep <- abs(tst$instRate - mu) <= outlier_sd * sig
    if (any(!keep, na.rm = TRUE)) tst <- tst[keep, , drop = FALSE]
  }

  # Protocol + percent change
  tst$protocol <- ifelse(tst$time <= 0, "baseline", "stimulus")
  bl_keep <- tst$time > bl_window[1] & tst$time <= bl_window[2]
  blMean <- mean(tst$instRate[bl_keep], na.rm = TRUE)
  tst$percent_change <- (tst$instRate / blMean) * 100

  # Plots
  if (isTRUE(plot_it)) {

    baseline_df <- spikeTTL$spike_puff_output
    wash_df     <- tst

    df_plot <- dplyr::bind_rows(
      baseline_df |> dplyr::select(time, percent_change) |> dplyr::mutate(condition = "baseline (spikeTTL)"),
      wash_df     |> dplyr::select(time, percent_change) |> dplyr::mutate(condition = "wash-in")
    )

    pred <- NA_character_
    if (!is.null(spikeTTL$selected_MD) && "assigned_subclass" %in% names(spikeTTL$selected_MD)) {
      pred <- spikeTTL$selected_MD$assigned_subclass[1]
    }

    ggplot2::ggplot(df_plot, ggplot2::aes(x = time, y = percent_change, colour = condition)) +
      ggplot2::geom_hline(yintercept = 100, linetype = 3) +
      ggplot2::geom_vline(xintercept = 0, linetype = 2) +
      ggplot2::geom_line(na.rm = TRUE) +
      ggplot2::labs(
        title = paste(cell, iM0$Exp, pred, sep = " - "),
        x = "time (s)",
        y = "% Change",
        linetype = NULL
      ) +
      ggplot2::theme_bw()


  }


  # Add metadata
  tst$cell_name <- iMD$cell_name
  tst$predicted_subclass <- iMD$predicted_subclass
  tst$assigned_subclass <- iMerged$assigned_subclass[!is.na(iMerged$assigned_subclass)][1]
  tst$Species <- iMD$Species
  tst$blockers <- iMD$Fast_blockers
  tst$bath <- iM0$Bath_Drug
  tst$puff <- iM0$Puff_drug
  tst$DFP <- iMD$Depth_from_pia
  tst$LIMs_depth <- iMD$LIMS_depth
  tst$expCon <- ifelse(is.na(iM0$Exp), "Standard_Puff", iM0$Exp)
  tst$response <- iM0$Effect

  # name sweep numbers
  names(sN) <- names(sweep_peaks)

  out <- list(
    cell = cell,
    md = iMD,
    selected_MD = iMerged,
    pro_testSpike = iM0,
    spike_TTL_sweep_numbers = sN,
    ttlTime = ttlTime,
    sweep_duration = dur,
    sweep_peaks = sweep_peaks,
    spike_puff_output = tst,
    blMean = blMean
  )

  if (isTRUE(write_csv)) {
    utils::write.csv(tst,
                     file.path(x$rd, paste0(cell, "-srtDrugWash.csv")),
                     row.names = FALSE)
  }

  out
}
