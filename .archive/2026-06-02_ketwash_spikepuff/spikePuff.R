spikePuff <- function(x,
                      mMD = NULL,
                      mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                      mMD_sheet = "Master_Macaque",
                      ttl_default = 5,
                      outlier_sd = 5,
                      bl_window = c(-5, 0),
                      write_csv = TRUE,
                      plot_it = TRUE,
                      ## NEW BEHAVIOR
                      rehydrate = TRUE,
                      save_srt = FALSE) {
  # -----------------------------
  # Load / accept srt
  # -----------------------------
  if (is.character(x)) {
    srt <- load5HT(x, tag = "-srt.rds")
  } else {
    srt <- x
  }
  cell <- srt$cell

  # -----------------------------
  # Load metadata tables (cheap)
  # -----------------------------
  MD <- headBCI::sheets$MD

  if (is.null(mMD)) {
    mMD <- readODS::read_ods(mMD_path, sheet = mMD_sheet)
    mMD <- mMD[!is.na(mMD$cell_name), , drop = FALSE]
  }

  # -----------------------------
  # Helpers (simple, editor-friendly)
  # -----------------------------
  pick_iM0 <- function(iM0) {
    if (nrow(iM0) > 1) {
      if ("Rank" %in% names(iM0)) {
        r1 <- which(iM0$Rank == 1 | iM0$Rank == "1")
        if (length(r1) > 0) {
          iM0 <- iM0[r1[1], , drop = FALSE]
        } else {
          bls <- which(grepl("blSpike", iM0$Rank))
          if (length(bls) > 0)
            iM0 <- iM0[bls[1], , drop = FALSE]
        }
      }
      if (nrow(iM0) > 1)
        iM0 <- iM0[1, , drop = FALSE]
    }
    iM0
  }

  refresh_md <- function() {
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

    m_idx <- which(grepl(cell, mMD$cell_name))
    if (length(m_idx) == 0) {
      message("No Master_UI_data match for cell: ", cell)
      return(NULL)
    }
    iMerged <- mMD[m_idx, , drop = FALSE]

    sp_idx <- which(grepl("spikeTTL", iMerged$stimulus_description))
    if (length(sp_idx) == 0) {
      message("No spikeTTL data for use")
      return(NULL)
    }
    iM0 <- iMerged[sp_idx, , drop = FALSE]
    iM0 <- pick_iM0(iM0)

    if (!("keep" %in% names(iM0))) {
      message("Missing 'keep' column in Master sheet subset.")
      return(NULL)
    }
    if (is.na(iM0$keep)) {
      message("check keep value. Keep is set NA. Need logical inference for analysis")
      return(NULL)
    }
    iM0 <- iM0[isTRUE(iM0$keep), , drop = FALSE]
    if (nrow(iM0) == 0) {
      message("No rows to keep for spikeTTL protocol in Master sheet subset.")
      return(NULL)
    }

    list(iMD = iMD,
         iMerged = iMerged,
         iM0 = iM0)
  }

  recompute_from_raw <- function(raw_trace, md_list) {
    iMD <- md_list$iMD
    iMerged <- md_list$iMerged
    iM0 <- md_list$iM0

    raw_trace <- raw_trace[!is.na(raw_trace$time) &
                             !is.na(raw_trace$instRate), , drop = FALSE]

    # outlier filter
    tst <- raw_trace
    mu  <- mean(tst$instRate, na.rm = TRUE)
    sig <- stats::sd(tst$instRate, na.rm = TRUE)
    if (is.finite(sig) && sig > 0) {
      keep <- abs(tst$instRate - mu) <= outlier_sd * sig
      if (any(!keep, na.rm = TRUE))
        tst <- tst[keep, , drop = FALSE]
    }

    # protocol
    tst$protocol <- ifelse(tst$time <= 0, "baseline", "stimulus")

    # baseline mean + percent change
    bl_keep <- tst$time > bl_window[1] & tst$time <= bl_window[2]
    blMean  <- mean(tst$instRate[bl_keep], na.rm = TRUE)
    tst$percent_change <- (tst$instRate / blMean) * 100

    # metadata columns (drop NA in assigned_subclass)
    assigned_sub <- iMerged$assigned_subclass
    assigned_sub <- assigned_sub[!is.na(assigned_sub)]
    assigned_sub <- if (length(assigned_sub) == 0)
      NA_character_
    else
      assigned_sub[1]

    tst$cell_name <- iMD$cell_name
    tst$predicted_subclass <- iMD$predicted_subclass
    tst$assigned_subclass <- assigned_sub
    tst$Species <- iMD$Species
    tst$blockers <- iMD$Fast_blockers
    tst$bath <- iM0$Bath_Drug
    tst$puff <- iM0$Puff_drug
    tst$DFP <- iMD$Depth_from_pia
    tst$LIMs_depth <- iMD$LIMS_depth
    tst$expCon <- ifelse(is.na(iM0$Exp), "Standard_Puff", iM0$Exp)
    tst$response <- iM0$Effect
    tst$keep <- isTRUE(iM0$keep)

    list(tst = tst, blMean = blMean)
  }

  # -----------------------------
  # Refresh metadata now
  # -----------------------------
  md_list <- refresh_md()
  if (is.null(md_list))
    return(NULL)

  iMD <- md_list$iMD
  iMerged <- md_list$iMerged
  iM0 <- md_list$iM0

  ttlTime <- ifelse(is.na(iM0$TTLtime), ttl_default, iM0$TTLtime)

  # -----------------------------
  # REHYDRATE (use srt$dfs$spikeTTL as cache)
  # -----------------------------
  has_cache <- !is.null(srt$dfs$spikeTTL) &&
    is.list(srt$dfs$spikeTTL) &&
    !is.null(srt$dfs$spikeTTL$spike_puff_output)

  if (isTRUE(rehydrate) && isTRUE(has_cache)) {
    # raw source preference: raw_trace (if you've cached it) else time/instRate from spike_puff_output
    raw_trace <- srt$dfs$spikeTTL$raw_trace
    if (is.null(raw_trace) ||
        !all(c("time", "instRate") %in% names(raw_trace))) {
      spo <- srt$dfs$spikeTTL$spike_puff_output
      raw_trace <- spo[, intersect(c("time", "instRate"), names(spo)), drop = FALSE]
      if (!all(c("time", "instRate") %in% names(raw_trace))) {
        message("Cannot rehydrate: cached spike_puff_output lacks time/instRate.")
        return(NULL)
      }
    }

    rec <- recompute_from_raw(raw_trace, md_list)
    tst <- rec$tst
    blMean <- rec$blMean

    # optional plot
    if (isTRUE(plot_it)) {
      assigned_all <- iMerged$assigned_subclass
      assigned_all <- assigned_all[!is.na(assigned_all)]
      main_txt <- paste(cell, paste(unique(assigned_all), collapse = ","), iM0$Puff_drug, sep = " - ")
      graphics::plot(
        tst$time,
        tst$percent_change,
        main = main_txt,
        xlab = "time",
        ylab = "% Change"
      )
    }

    # overwrite cached outputs + refreshed metadata
    srt$dfs$spikeTTL$cell <- cell
    srt$dfs$spikeTTL$md <- iMD
    srt$dfs$spikeTTL$selected_MD <- iMerged
    srt$dfs$spikeTTL$pro_spikeTTL <- iM0
    srt$dfs$spikeTTL$ttlTime <- ttlTime
    srt$dfs$spikeTTL$blMean <- blMean
    srt$dfs$spikeTTL$spike_puff_output <- tst
    srt$dfs$spikeTTL$raw_trace <- raw_trace

    if (isTRUE(write_csv)) {
      utils::write.csv(tst, file.path(srt$rd, paste0(cell, "-srtPuff.csv")), row.names = FALSE)
    }

    if (isTRUE(save_srt)) {
      # If your project has a canonical save helper, swap this line.
      saveRDS(srt, file.path(srt$rd, paste0(cell, "-srt.rds")))
    }

    return(srt$dfs$spikeTTL)
  }

  # -----------------------------
  # HEAVY BUILD (compute then store in srt$dfs$spikeTTL)
  # -----------------------------

  # Parse sweeps
  pro <- as.character(iM0$stimulus_description)
  spike_puff_sweeps <- as.character(iM0$wrapped_sweeps)
  spike_puff_sweeps <- gsub("\\[|\\]", "", spike_puff_sweeps)

  sp <- strsplit(spike_puff_sweeps, ",")[[1]]
  sp <- trimws(sp)
  sN <- as.numeric(sp) + 1

  # heavy: extract NWB sweeps
  nwb <- nphys::extractNWB(iMD$nwb_file, sweeps = sN)

  if (is.null(names(nwb$sweeps))) {
    nwb$sweeps <- as.data.frame(nwb$sweeps)
  }
  sweep_names <- names(nwb$sweeps)

  # durations
  dur <- vapply(sweep_names, function(nm) {
    sw <- nwb$sweeps[[nm]]
    sw <- sw[!is.na(sw)]
    sr <- nwb$sampling_rate[[nm]]
    tail(nphys::convertIndex(sw, sr), 1)
  }, numeric(1))

  # heavy: APanalysis
  sweep_peaks <- lapply(sweep_names, function(nm) {
    sweep <- nwb$sweeps[[nm]]
    sweep <- sweep[!is.na(sweep)]
    sr <- nwb$sampling_rate[[nm]]
    nphys::APanalysis(sweep, rate = sr)
  })
  names(sweep_peaks) <- sweep_names

  # name sweeps by protocol chunks
  pro_parts <- unlist(strsplit(pro, "\\+"))
  if (length(pro_parts) != length(sweep_peaks)) {
    names(sweep_peaks) <- c(pro_parts, rep("trail", max(
      0, length(sweep_peaks) - length(pro_parts)
    )))
  } else {
    names(sweep_peaks) <- pro_parts
  }
  names(dur) <- names(sweep_peaks)

  # Build aligned data frames
  dfBl <- dfTTL <- dfTrail <- NULL

  # baseline
  if ("baseline" %in% names(sweep_peaks)) {
    fp <- as.data.frame(sweep_peaks[[which(names(sweep_peaks) == "baseline")[1]]])
    fp <- fp[-nrow(fp), , drop = FALSE]
    maxpk <- as.numeric(dur[which(names(dur) == "baseline")[1]])
    dfBl <- data.frame(time = fp$epoch_t - (maxpk + ttlTime),
                       instRate = fp$instRate)
  }

  # spikeTTL
  if ("spikeTTL" %in% names(sweep_peaks)) {
    fp <- as.data.frame(sweep_peaks[[which(names(sweep_peaks) == "spikeTTL")[1]]])
    fp <- fp[-nrow(fp), , drop = FALSE]
    dfTTL <- data.frame(time = fp$epoch_t - ttlTime,
                        instRate = fp$instRate)
  } else {
    message("No spikeTTL sweep_peaks found after naming; returning NULL.")
    return(NULL)
  }

  # trail
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

      fp_all <- do.call(rbind, fps)
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

  # Combine into raw trace (store for future rehydrate)
  tst_raw <- dfTTL
  if (!is.null(dfBl))
    tst_raw <- rbind(dfBl, tst_raw)
  if (!is.null(dfTrail))
    tst_raw <- rbind(tst_raw, dfTrail)

  # Cheap recompute
  rec <- recompute_from_raw(tst_raw, md_list)
  tst <- rec$tst
  blMean <- rec$blMean

  # plot
  if (isTRUE(plot_it)) {
    assigned_all <- iMerged$assigned_subclass
    assigned_all <- assigned_all[!is.na(assigned_all)]
    main_txt <- paste(cell, paste(unique(assigned_all), collapse = ","), iM0$Puff_drug, sep = " - ")
    graphics::plot(
      tst$time,
      tst$percent_change,
      main = main_txt,
      xlab = "time",
      ylab = "% Change"
    )
  }

  # rename sweep numbers/names for output
  names(sN) <- names(sweep_peaks)
  names(sweep_names) <- names(sweep_peaks)

  out <- list(
    cell = cell,
    md = iMD,
    selected_MD = iMerged,
    pro_spikeTTL = iM0,
    sweep_names = sweep_names,
    spike_TTL_sweep_numbers = sN,
    ttlTime = ttlTime,
    sweep_duration = dur,
    sweep_peaks = sweep_peaks,
    raw_trace = tst_raw,
    # <--- key for rehydrate
    spike_puff_output = tst,
    blMean = blMean
  )

  if (isTRUE(write_csv)) {
    utils::write.csv(tst, file.path(srt$rd, paste0(cell, "-srtPuff.csv")), row.names = FALSE)
  }

  if (isTRUE(save_srt)) {
    # store in srt cache
    srt$dfs$spikeTTL <- out
    saveRDS(srt, file.path(srt$rd, paste0(cell, "-srt.rds")))
  }

  out
}
