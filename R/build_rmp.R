#' Build RMP protocol structures from tpRMP sweeps
#'
#' Drop-in analogue of \code{build_asp()}, but for tpRMP (membrane potential).
#'
#' Key behaviors:
#' - selects protocol rows where stimulus_description contains "tpRMP_puff"
#' - chunk mapping uses "tpRMPbl" as baseline and preserves "trail"
#' - stores outputs in \code{srt$dfs$rmp}
#'
#' @param x A \code{srt} object (list) or a character cell ID (will be loaded via \code{load5HT()}).
#' @param mMD Optional Master sheet data.frame; if NULL it is read from \code{mMD_path}.
#' @param mMD_path Path to Master_UI_data.ods
#' @param mMD_sheet Sheet name inside the ODS
#' @param puff_default Default puff time if TTLtime is missing.
#' @param outlier_sd Outlier filter threshold in SDs (tpRMP).
#' @param bl_window Baseline window on aligned time axis, default c(-5, 0).
#' @param save_srt Save updated srt back to disk.
#' @param write_csv Write aligned time-series CSV per protocol.
#' @param store_level "light" or "full" (controls plot artifact retention).
#' @param csv_tag Suffix for CSV outputs.
#'
#' @return The \code{rmp} container (\code{srt$dfs$rmp}).
#' @export
build_rmp <- function(x,
                      mMD = NULL,
                      mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                      mMD_sheet = "Master_Macaque",
                      puff_default = 5,
                      outlier_sd = 5,
                      bl_window = c(-5, 0),
                      save_srt = TRUE,
                      write_csv = FALSE,
                      store_level = c("light", "full"),
                      csv_tag = "-rmp.csv") {

  store_level <- match.arg(store_level)

  # ---- Load / accept srt ----
  if (is.character(x)) {
    srt <- load5HT(x, tag = "-srt.rds")
  } else {
    srt <- x
  }
  cell <- srt$cell

  # ---- Load metadata ----
  MD <- headBCI::sheets$MD

  if (is.null(mMD)) {
    mMD <- readODS::read_ods(mMD_path, sheet = mMD_sheet)
    mMD <- mMD[!is.na(mMD$cell_name), , drop = FALSE]
  }

  md_idx <- which(as.character(MD$cell_name) == cell)
  if (length(md_idx) == 0) {
    message("No MD match for cell: ", cell)
    return(NULL)
  }
  iMD <- MD[md_idx[1], , drop = FALSE]

  m_idx <- which(as.character(mMD$cell_name) == cell)
  if (length(m_idx) == 0) {
    message("No Master_UI_data match for cell: ", cell)
    return(NULL)
  }
  iMerged <- mMD[m_idx, , drop = FALSE]

  # per-cell fills (match your ASP behavior)
  iMerged$assigned_subclass <- fill_first_nonempty(iMerged$assigned_subclass)
  iMerged$TTLtime <- fill_locf(iMerged$TTLtime, default = puff_default, type = "numeric")
  iMerged$Exp[is.na(iMerged$Exp)] <- "standard_puff"

  # ---- Identify RMP protocols ----
  stim <- as.character(iMerged$stimulus_description)
  idx <- which(grepl("tpRMP_puff", stim, fixed = TRUE))
  if (length(idx) == 0) {
    message("No tpRMP_puff protocols found for cell: ", cell)
    return(NULL)
  }

  iP <- iMerged[idx, , drop = FALSE]
  if (!("keep" %in% names(iP))) {
    message("Missing 'keep' column in Master sheet subset.")
    return(NULL)
  }

  keep_vec <- as_keep_logical(iP$keep)
  iP <- iP[keep_vec, , drop = FALSE]
  if (nrow(iP) == 0) {
    message("No tpRMP_puff rows with keep == TRUE for cell: ", cell)
    return(NULL)
  }

  # ---- init container ----
  if (is.null(srt$dfs$rmp) || !is.list(srt$dfs$rmp)) srt$dfs$rmp <- list()
  srt$dfs$rmp$by_protocol <- list()
  smry_rows <- list()

  # ---- loop protocols ----
  #k = 1

  for (k in seq_len(nrow(iP))) {

    iM0 <- iP[k, , drop = FALSE]
    puffTime <- ifelse(is.na(iM0$TTLtime), puff_default, iM0$TTLtime)
    sN <- parse_wrapped_sweeps(iM0$wrapped_sweeps)

    nwb <- nphys::extractNWB(iMD$nwb_file, sweeps = sN)
    if (is.null(names(nwb$sweeps))) nwb$sweeps <- as.data.frame(nwb$sweeps)
    sweep_names <- names(nwb$sweeps)

    dur <- vapply(sweep_names, function(nm) {
      sw <- nwb$sweeps[[nm]]
      sw <- sw[!is.na(sw)]
      sr <- nwb$sampling_rate[[nm]]
      tail(nphys::convertIndex(sw, sr), 1)
    }, numeric(1))

    # ---- build sweep_objs with epochDF fallback ----
    prev_epDF <- NULL
    sweep_objs <- vector("list", length(sweep_names))
    names(sweep_objs) <- sweep_names
    #ii = 1
    for (ii in seq_along(sweep_names)) {

      nm <- sweep_names[ii]
      sw <- nwb$sweeps[[nm]]
      sr <- nwb$sampling_rate[[nm]]
      df <- sweep_to_tpRMP_df(sw, sr)

      ep_source <- NA_character_
      epDF <- NULL

      # 1) epochDF from comments
      ei <- paste0(nm, ".comment")
      if (!is.null(nwb$epochs) && !is.null(nwb$epochs[[ei]]) && "epochDF" %in% names(nwb$epochs[[ei]])) {
        epDF <- tryCatch(as.data.frame(nwb$epochs[[ei]][["epochDF"]]), error = function(e) NULL)
      }
      if (is_epochDF_usable(epDF, "T0", "T1")) {
        ep_source <- "epochDF"
      } else {
        epDF <- NULL
      }

      # 2) fallback prev
      if (is.null(epDF) && is_epochDF_usable(prev_epDF, "T0", "T1")) {
        epDF <- prev_epDF
        ep_source <- "prev_epochDF"
      }

      # 3) fallback stimulus
      if (is.null(epDF)) {
        stim_vec <- get_stim_for_sweep(nwb, nm)
        #plot(stim_vec, type = "l")
        if (is.null(stim_vec)) {
          df$pulse_state <- "pulse_off"
          sweep_objs[[nm]] <- list(df = df, seg = data.frame(), sr = sr, ep_source = "no_epochs_no_stim")
          next
        }

        epDF <- make_epochDF_from_stimulus(stim_vec, sr = sr, pulse_polarity = "auto")
        if (!is_epochDF_usable(epDF, "T0", "T1")) {
          df$pulse_state <- "pulse_off"
          sweep_objs[[nm]] <- list(df = df, seg = data.frame(), sr = sr, ep_source = "stim_no_pulses")
          next
        }

        ep_source <- "stimulus"
      }

      epDF_clean <- clean_epochDF(epDF, "T0", "T1")
      if (is.null(epDF_clean) || nrow(epDF_clean) == 0) {
        df$pulse_state <- "pulse_off"
        sweep_objs[[nm]] <- list(df = df, seg = data.frame(), sr = sr, ep_source = "ep_clean_empty")
        next
      }

      seg <- tryCatch(
        order_epoch_segments(epDF_clean, sampling_rate = sr, t0_col = "T0", t1_col = "T1"),
        error = function(e) data.frame()
      )

      if (!is.data.frame(seg) || nrow(seg) == 0) {
        df$pulse_state <- "pulse_off"
        sweep_objs[[nm]] <- list(df = df, seg = data.frame(), sr = sr, ep_source = "seg_empty")
        next
      }

      df <- add_pulse_state_time(df, seg,
                                 df_time_col = "epoch_t",
                                 seg_t0_col = "T0",
                                 seg_t1_col = "T1",
                                 seg_subtype_col = "subtype")

      prev_epDF <- epDF
      sweep_objs[[nm]] <- list(df = df, seg = seg, sr = sr, ep_source = ep_source)
    }

    # ---- chunk parse + build raw ----
    chunk_info  <- parse_stim_chunks_rmp(iM0$stimulus_description, length(sweep_objs))
    chunk_names <- chunk_info$parts

    built <- build_raw_trace_rmp(sweep_objs, dur, chunk_names, puffTime)
    raw_trace <- built$raw_trace
    segment_bounds <- built$segment_bounds

    protocol_key <- make_protocol_key_rmp(iM0, iMD, chunk_info, length(sweep_objs), puff_default = puff_default)

    if (is.null(raw_trace) || nrow(raw_trace) == 0) {
      message("RMP: empty raw trace for ", cell, " | ", protocol_key)
      next
    }

    rec <- recompute_from_raw_rmp(raw_trace, outlier_sd = outlier_sd, bl_window = bl_window)
    tst <- rec$tst
    blMean <- rec$blMean

    # ---- analysis ----
    proto_analysis <- analyze_rmp_protocol_set(
      tst,
      downsample_for_plot = (store_level == "full"),
      make_plots = TRUE,
      return_level = ifelse(store_level == "full", "full", "light")
    )



    protocol_map <- list(
      chunk_names = chunk_names,
      baseline_pos = chunk_info$baseline_pos,
      trail_pos = chunk_info$trail_pos,
      puff_sweep_positions = chunk_info$puff_pos,
      n_puff = length(chunk_info$puff_pos)
    )

    out_k <- list(
      protocol_key = protocol_key,
      cell = cell,
      md = iMD,
      selected_MD = iMerged,
      pro_row = iM0,
      protocol_map = protocol_map,
      puffTime = puffTime,
      sweep_duration = dur,
      segment_bounds = segment_bounds,
      blMean = blMean,
      epoch_sources = table(vapply(sweep_objs, function(x) x$ep_source %||% NA_character_, character(1)))
    )

    out_k$rmp_protocol_analysis <- list(
      chunk_summary = proto_analysis$chunk_summary,
      bl_mean = proto_analysis$bl_mean,
      trail_mean = proto_analysis$trail_mean,
      pulse = proto_analysis$pulse
    )

    if (store_level == "full") {
      out_k$rmp_protocol_analysis$plot_trace_df <- proto_analysis$plot_trace_df
      out_k$plots <- list(
        rmp_trace = proto_analysis$plot_trace,
        pulse_features = if (!is.null(proto_analysis$pulse)) proto_analysis$pulse$plot_features else NULL
      )
    } else {
      out_k$plots <- NULL
    }

    srt$dfs$rmp$by_protocol[[protocol_key]] <- out_k

    smry_rows[[protocol_key]] <- data.frame(
      protocol_key = protocol_key,
      stimulus_description = trim1(iM0$stimulus_description),
      n_sweeps = length(sweep_objs),
      n_puff = protocol_map$n_puff,
      puff_positions = paste(protocol_map$puff_sweep_positions, collapse = ","),
      wrapped_sweeps = trim1(iM0$wrapped_sweeps),
      puffTime = puffTime,
      Bath_Drug = safe_chr(iMD$bath_drug),
      Puff_drug = safe_chr(iMD$puff_drug),
      Exp = safe_chr(iM0$Exp),
      Effect = safe_chr(iM0$Effect),
      stringsAsFactors = FALSE
    )

    if (isTRUE(write_csv)) {
      fn <- file.path(srt$rd, paste0(cell, "_", protocol_key, csv_tag))
      utils::write.csv(tst, fn, row.names = FALSE)
    }
  }

  # ---- header fields ----
  srt$dfs$rmp$cell <- cell
  srt$dfs$rmp$md <- iMD
  srt$dfs$rmp$selected_MD <- iMerged
  srt$dfs$rmp$protocol_keys <- names(srt$dfs$rmp$by_protocol)
  srt$dfs$rmp$smry <- if (length(smry_rows)) do.call(rbind, smry_rows) else NULL

  if (isTRUE(save_srt)) {
    saveRDS(srt, file.path(srt$rd, paste0(cell, "-srt.rds")))
  }

  srt$dfs$rmp
}
