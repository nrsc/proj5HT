#' Build ASP (aligned spike protocol) structures from spikeTTL / TTLstack sweeps
#'
#' Constructs an ASP container from a per-cell `srt` object (or a `cell_name` that can
#' be loaded via `load5HT()`), identifying protocols in the Master UI sheet whose
#' `stimulus_description` indicates `spikeTTL` or `TTLstack` and whose `keep` flag is
#' truthy. For each kept protocol row, the function:
#' \enumerate{
#'   \item extracts the listed sweeps from the cell NWB file,
#'   \item computes per-sweep spike-rate traces via `nphys::APanalysis()`,
#'   \item builds an aligned local-time trace (`raw_trace`) where TTL-containing sweeps
#'         are shifted so TTL onset lands near 0 seconds,
#'   \item recomputes baseline mean and percent-of-baseline metrics within `bl_window`,
#'   \item stores per-protocol outputs in `asp$by_protocol` under a deterministic
#'         `protocol_key`.
#' }
#'
#' The returned object is designed for downstream plotting (e.g., concatenated protocol
#' views, per-protocol trace + heat strip) and includes `segment_bounds` with sweep-level
#' start/end times and (when possible) a computed `ttl_onset_time`.
#'
#' @param x Either a character `cell_name` (loaded via `load5HT(cell, tag="-srt.rds")`)
#'   or an existing `srt` nnested list that contains at least `srt$cell`, `srt$rd`,
#'   and `srt$dfs`.
#' @param mMD Optional data.frame containing the Master UI metadata (must include at
#'   least `cell_name`, `stimulus_description`, `wrapped_sweeps`, `keep`, and typically
#'   columns like `Bath_Drug`, `Puff_drug`, `TTLtime`, `Exp`, `Effect`,
#'   `assigned_subclass`). If `NULL`, it is read from `mMD_path` / `mMD_sheet`.
#' @param mMD_path Character. Path to the Master UI `.ods` file used when `mMD` is
#'   not supplied.
#' @param mMD_sheet Character. Sheet name within the `.ods` file used when `mMD` is
#'   not supplied.
#' @param ttl_default Numeric. Default TTL onset time (seconds) used when the protocol
#'   row lacks a valid `TTLtime`.
#' @param outlier_sd Numeric. Outlier threshold (in SD units) applied to `instRate`
#'   within each protocol trace before recomputing baseline/percent change.
#'   Points where `|instRate - mean(instRate)| > outlier_sd * sd(instRate)` are dropped.
#' @param bl_window Numeric length-2 vector. Baseline window in seconds (e.g. `c(-5,0)`)
#'   used to compute `blMean` and `percent_change`.
#' @param save_srt Logical. If `TRUE`, saves the updated `srt` object back to disk as
#'   `file.path(srt$rd, paste0(cell, "-srt.rds"))`.
#' @param plot_it Logical. If `TRUE`, produces a quick base R plot of `percent_change`
#'   vs time for each protocol during building (primarily for debugging).
#' @param write_csv Logical. If `TRUE`, writes a per-protocol CSV of `spike_puff_output`
#'   to `srt$rd` using `csv_tag`.
#' @param csv_tag Character. Filename suffix appended to per-protocol CSV exports when
#'   `write_csv = TRUE` (default `"-asp.csv"`).
#'
#' @return An ASP list (also stored at `srt$dfs$asp`) with components:
#' \describe{
#'   \item{`cell`}{Cell name / identifier.}
#'   \item{`md`}{Row from `headBCI::sheets$MD` corresponding to `cell`.}
#'   \item{`selected_MD`}{Subset of Master UI rows (`mMD`) corresponding to `cell`
#'     (may include multiple rows).}
#'   \item{`protocol_keys`}{Character vector of protocol keys created.}
#'   \item{`smry`}{Data.frame summary table of protocols (n sweeps, n TTL sweeps, sweep
#'     positions, wrapped_sweeps, TTLtime, Bath/Puff/Exp/Effect, etc.).}
#'   \item{`by_protocol`}{Named list of per-protocol objects keyed by `protocol_key`. Each
#'     protocol entry contains:
#'     \describe{
#'       \item{`protocol_key`}{Deterministic key describing protocol contents/metadata.}
#'       \item{`pro_row`}{The Master UI row used to build this protocol.}
#'       \item{`protocol_map`}{Chunk labels and positional indices (baseline, trail, TTL sweeps).}
#'       \item{`ttlTime`}{TTL time (seconds) used for TTL sweep alignment.}
#'       \item{`sweep_duration`}{Per-sweep durations (seconds).}
#'       \item{`sweep_peaks`}{List of `nphys::APanalysis()` outputs per sweep.}
#'       \item{`segment_bounds`}{Sweep-level bounds with `start_time`, `end_time`, and
#'         `ttl_onset_time` where computable.}
#'       \item{`raw_trace`}{Aligned local-time trace (`time`, `instRate`) built by stitching sweeps.}
#'       \item{`spike_puff_output`}{Processed per-timepoint output including `protocol` flag,
#'         `percent_change`, and metadata columns (cell, subclass, species, blockers, bath/puff, etc.).}
#'       \item{`blMean`}{Baseline mean instRate used to compute `percent_change`.}
#'     }}
#' }
#'
#' @details
#' **Protocol identification:** protocols are selected from the Master UI sheet by
#' matching `stimulus_description` containing `spikeTTL` or `TTLstack`, then filtering
#' rows by a permissive truthy parser of the `keep` column (accepts logicals, "yes"/"x",
#' "keep", non-zero numbers, etc.).
#'
#' **Sweep parsing:** `wrapped_sweeps` is parsed into numeric sweep ids and incremented
#' by 1 (to match NWB sweep indexing used here).
#'
#' **Chunk mapping:** `stimulus_description` is split on `"+"` to produce sweep-level
#' chunk names; `TTLstack` chunks are expanded to unique labels (`TTLstack_1`, `TTLstack_2`, ...).
#' TTL-containing sweeps are those labeled `spikeTTL` or `TTLstack_*`.
#'
#' **Time alignment:** the first TTL-containing sweep defines the zero reference for
#' concatenation; TTL sweeps are additionally shifted by `ttlTime` so TTL onset occurs
#' near 0 seconds within those sweeps. `segment_bounds$ttl_onset_time` is computed by
#' snapping to the sample closest to `ttlTime` within each TTL sweep.
#'
#' **Recomputation:** after assembling the raw trace, `instRate` is outlier-filtered,
#' `blMean` is computed over `bl_window`, and `percent_change` is computed as
#' `(instRate / blMean) * 100` when possible (otherwise `NA`).
#'
#' @examples
#' \dontrun{
#' # Build from a cell name (loads -srt.rds), using default Master sheet
#' asp <- build_asp("Q21.26.014.1A.21.02")
#'
#' # Build using a preloaded srt object and a preloaded Master sheet
#' mMD <- readODS::read_ods("~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
#'                         sheet = "Master_Macaque")
#' asp <- build_asp(srt, mMD = mMD, save_srt = FALSE)
#'
#' # Debug plots and export per-protocol CSVs
#' asp <- build_asp("Q21.26.014.1A.21.02", plot_it = TRUE, write_csv = TRUE)
#'
#' x = srt
#' }
#'
#' @export
build_asp <- function(x,
                      mMD = NULL,
                      mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                      mMD_sheet = "Master_Macaque",
                      ttl_default = 5,
                      outlier_sd = 5,
                      bl_window = c(-5, 0),
                      save_srt = TRUE,
                      plot_it = FALSE,
                      write_csv = TRUE,
                      csv_tag = "-asp.csv") {



  # -----------------------------
  # Helper functions
  # -----------------------------
  trim1 <- function(z) trimws(gsub("\\s+", " ", as.character(z)))

  safe_chr <- function(z) ifelse(is.na(z), "", as.character(z))

  # ============================================
  # Helpers for filling Master UI columns
  # ============================================

  # Repeat the first non-empty (non-NA, non-"") value for the whole column
  # (good for per-cell constants like assigned_subclass)
  fill_first_nonempty <- function(x, default = NA_character_) {
    x_chr <- trimws(as.character(x))
    x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "nan")] <- NA_character_

    use <- x_chr[!is.na(x_chr)]
    use <- if (length(use) > 0) use[1] else as.character(default)

    rep(use, length(x_chr))
  }

  # Carry last observation forward (LOCF) with a default for leading missing values.
  # - If the first value is missing/empty => set to default
  # - Then carry-forward through the vector
  # - Remaining NAs (if any) get default
  fill_locf <- function(x, default, type = c("character", "numeric")) {
    type <- match.arg(type)

    x_chr <- trimws(as.character(x))
    x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "nan")] <- NA_character_

    if (type == "numeric") {
      x_out <- suppressWarnings(as.numeric(x_chr))
      default_val <- as.numeric(default)
    } else {
      x_out <- x_chr
      default_val <- as.character(default)
    }

    if (length(x_out) == 0L) return(x_out)

    # leading NA -> default
    if (is.na(x_out[1])) x_out[1] <- default_val

    # carry-forward
    for (ii in seq_along(x_out)) {
      if (ii == 1) next
      if (is.na(x_out[ii])) x_out[ii] <- x_out[ii - 1]
    }

    # any remaining NA -> default
    x_out[is.na(x_out)] <- default_val
    x_out
  }

  as_keep_logical <- function(x) {
    if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))

    x0 <- tolower(trimws(as.character(x)))
    x0[x0 %in% c("", "na", "nan")] <- NA_character_

    out <- rep(FALSE, length(x0))

    # TRUE-like
    out[!is.na(x0) & x0 %in% c("true","t","1","yes","y","keep","k","x")] <- TRUE
    # FALSE-like
    out[!is.na(x0) & x0 %in% c("false","f","0","no","n","drop","d")] <- FALSE

    # numeric fallback
    suppressWarnings({
      xn <- as.numeric(x0)
      out[!is.na(xn)] <- xn[!is.na(xn)] != 0
    })

    out[is.na(x0)] <- FALSE
    out
  }

  parse_wrapped_sweeps <- function(w) {
    w <- gsub("\\[|\\]", "", as.character(w))
    sp <- strsplit(w, ",")[[1]]
    sp <- trimws(sp)
    as.numeric(sp) + 1
  }

  # Define chunk names + TTL sweep positions from stim_desc
  parse_stim_chunks <- function(stim_desc, n_sweeps) {
    parts0 <- unlist(strsplit(as.character(stim_desc), "\\+"))
    parts0 <- trimws(parts0)
    parts0 <- parts0[parts0 != ""]

    if (length(parts0) != n_sweeps) {
      # conservative fallback
      if (length(parts0) == 0) parts0 <- rep("unknown", n_sweeps)
      if (length(parts0) < n_sweeps) parts0 <- c(parts0, rep("trail", n_sweeps - length(parts0)))
      if (length(parts0) > n_sweeps) parts0 <- parts0[seq_len(n_sweeps)]
    }

    # Expand TTLstack into unique sweep-level labels
    ttlstack_idx <- which(grepl("TTLstack", parts0, fixed = TRUE))
    parts <- parts0
    if (length(ttlstack_idx) > 0) {
      for (j in seq_along(ttlstack_idx)) {
        parts[ttlstack_idx[j]] <- paste0("TTLstack_", j)
      }
    }

    ttl_pos <- which(parts == "spikeTTL" | grepl("^TTLstack_", parts))
    list(
      parts = parts,
      ttl_pos = ttl_pos,
      baseline_pos = which(parts == "baseline"),
      trail_pos = which(parts == "trail")
    )
  }

  # Ensure schema: always have instRate & percent_change when possible
  ensure_instRate_percent_change <- function(tst, blMean = NULL) {
    if (is.null(tst) || !is.data.frame(tst)) return(tst)

    # instRate: if missing but percent_change + blMean exist, reconstruct
    if (!("instRate" %in% names(tst))) {
      if ("percent_change" %in% names(tst) &&
          !is.null(blMean) && is.finite(blMean) && blMean > 0) {
        tst$instRate <- (tst$percent_change / 100) * blMean
      } else {
        tst$instRate <- NA_real_
      }
    }

    # percent_change: if missing but instRate + blMean exist, compute
    if (!("percent_change" %in% names(tst))) {
      if ("instRate" %in% names(tst) &&
          !is.null(blMean) && is.finite(blMean) && blMean > 0) {
        tst$percent_change <- (tst$instRate / blMean) * 100
      } else {
        tst$percent_change <- NA_real_
      }
    }

    tst
  }

  recompute_from_raw <- function(raw_trace, iMD, iMerged, iM0) {
    raw_trace <- raw_trace[!is.na(raw_trace$time) & !is.na(raw_trace$instRate), , drop = FALSE]

    # outlier filter
    tst <- raw_trace
    mu  <- mean(tst$instRate, na.rm = TRUE)
    sig <- stats::sd(tst$instRate, na.rm = TRUE)
    if (is.finite(sig) && sig > 0) {
      keep <- abs(tst$instRate - mu) <= outlier_sd * sig
      if (any(!keep, na.rm = TRUE)) tst <- tst[keep, , drop = FALSE]
    }

    # protocol flag
    tst$protocol <- ifelse(tst$time <= 0, "baseline", "stimulus")

    # baseline mean + percent change
    bl_keep <- tst$time > bl_window[1] & tst$time <= bl_window[2]
    blMean  <- mean(tst$instRate[bl_keep], na.rm = TRUE)

    if (!is.finite(blMean) || blMean <= 0) {
      blMean <- NA_real_
      tst$percent_change <- NA_real_
    } else {
      tst$percent_change <- (tst$instRate / blMean) * 100
    }

    assigned_sub <- iMerged$assigned_subclass
    assigned_sub <- assigned_sub[!is.na(assigned_sub)]
    assigned_sub <- if (length(assigned_sub) == 0) NA_character_ else as.character(assigned_sub[1])


    # metadata columns
    tst$cell_name <- iMD$cell_name
    tst$predicted_subclass <- iMD$predicted_subclass
    tst$Species <- iMD$Species
    tst$blockers <- iMD$Fast_blockers
    tst$bath <- iMD$bath_drug
    tst$puff <- iMD$puff_drug
    tst$DFP <- iMD$Depth_from_pia
    tst$LIMs_depth <- iMD$LIMS_depth
    tst$assigned_subclass <- rep(assigned_sub, nrow(tst))
    tst$expCon <- safe_chr(as.character(iM0$Exp[1]))
    tst$response <- iM0$Effect
    #tst$keep <- isTRUE(iM0$keep)
    tst$keep <- as_keep_logical(iM0$keep)[1]


    # enforce schema (in case baseline failed)
    tst <- ensure_instRate_percent_change(tst, blMean = blMean)

    list(tst = tst, blMean = blMean)
  }

  make_protocol_key <- function(iM0, iMerged, chunk_info, n_sweeps) {

    stim_clean <- trim1(iM0$stimulus_description[1])
    typ <- if (grepl("TTLstack", stim_clean)) "TTLstack" else if (grepl("spikeTTL", stim_clean)) "spikeTTL" else "other"
    ttl <- ifelse(is.na(iM0$TTLtime[1]), ttl_default, iM0$TTLtime[1])

    key <- paste(
      typ,
      paste0("subc=", iM0$assigned_subclass[1]),
      paste0("nSweep=", n_sweeps),
      paste0("nTTL=", length(chunk_info$ttl_pos)),
      paste0("sweeps=", trim1(iM0$wrapped_sweeps[1])),
      paste0("ttlTime=", ttl),
      paste0("bath=", safe_chr(iMD$bath_drug[1])),
      paste0("puffdrug=", safe_chr(iMD$puff_drug[1])),
      paste0("exp=", safe_chr(iM0$Exp[1])),
      sep = "|"
    )

    gsub("[^A-Za-z0-9|=_\\-]+", "_", key)
  }

  # Build aligned raw trace + segment_bounds from sweep_peaks/dur/chunks
  build_raw_trace <- function(sweep_peaks, dur, chunk_names, ttlTime) {
    # name by chunk
    names(sweep_peaks) <- chunk_names
    names(dur) <- chunk_names

    peaks_to_df <- function(pk) {
      d <- as.data.frame(pk)
      if (nrow(d) >= 2) d <- d[-nrow(d), , drop = FALSE]
      data.frame(epoch_t = d$epoch_t, instRate = d$instRate)
    }

    dfs <- lapply(sweep_peaks, peaks_to_df)

    # continuous shifts by duration
    cum_shift <- c(0, cumsum(as.numeric(dur)))[seq_along(dfs)]

    # first TTL sweep defines "zero"
    first_ttl_idx <- which(chunk_names == "spikeTTL" | grepl("^TTLstack_", chunk_names))[1]
    if (is.na(first_ttl_idx)) first_ttl_idx <- 1L
    zero_offset <- cum_shift[first_ttl_idx]

    segment_bounds <- data.frame(
      sweep = seq_along(dfs),
      chunk = chunk_names,
      start_time = NA_real_,
      end_time = NA_real_,
      ttl_onset_time = NA_real_,   # NEW
      stringsAsFactors = FALSE
    )

    out_list <- vector("list", length(dfs))

    for (i in seq_along(dfs)) {
      d <- dfs[[i]]

      t_pre <- d$epoch_t + cum_shift[i] - zero_offset
      #t <- t_pre
      t <- t_pre - ttlTime

      out_list[[i]] <- data.frame(time = t, instRate = d$instRate)

      segment_bounds$start_time[i] <- min(t, na.rm = TRUE)
      segment_bounds$end_time[i]   <- max(t, na.rm = TRUE)
    }


    raw <- do.call(rbind, out_list)
    list(raw_trace = raw, segment_bounds = segment_bounds)
  }
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
  # Load metadata
  # -----------------------------
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
  if (length(md_idx) > 1) {
    message("Multiple MD matches for cell ", cell, " (using first).")
    md_idx <- md_idx[1]
  }
  iMD <- MD[md_idx, , drop = FALSE]


  m_idx <- which(as.character(mMD$cell_name) == cell)
  if (length(m_idx) == 0) {
    message("No Master_UI_data match for cell: ", cell)
    return(NULL)
  }
  if (length(m_idx) > 1) {
    message("Multiple Master_UI_data matches for cell ", cell, " (keeping all).")
  }
  iMerged <- mMD[m_idx, , drop = FALSE]

  # assigned_subclass: per-cell constant
  iMerged$assigned_subclass <- fill_first_nonempty(iMerged$assigned_subclass)

  # TTLtime: carry-forward, leading NA means 5
  iMerged$TTLtime <- fill_locf(iMerged$TTLtime, default = 5, type = "numeric")

  # Exp: your current rule (only replace NA with standard_puff)
  iMerged$Exp[is.na(iMerged$Exp)] <- "standard_puff"

  # If Exp ever needs the same "NA NA NA  Ket-wash-in Ket-wash-in" behavior, do this instead:
  # iMerged$Exp <- fill_locf(iMerged$Exp, default = "standard_puff", type = "character")



  # -----------------------------
  # Identify asp protocols: spikeTTL or TTLstack
  # -----------------------------
  stim <- as.character(iMerged$stimulus_description)
  idx <- which(grepl("spikeTTL", stim) | grepl("TTLstack", stim))
  if (length(idx) == 0) {
    message("No spikeTTL/TTLstack protocols found for cell: ", cell)
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
    message("No spikeTTL/TTLstack rows with keep == TRUE for cell: ", cell)
    return(NULL)
  }


  # -----------------------------
  # Init asp container
  # -----------------------------
  if (is.null(srt$dfs$asp) || !is.list(srt$dfs$asp)) srt$dfs$asp <- list()
  srt$dfs$asp$by_protocol <- list()
  smry_rows <- list()

  # -----------------------------
  # Loop all protocols
  # -----------------------------
  for (k in seq_len(nrow(iP))) {

    iM0 <- iP[k, , drop = FALSE]
    ttlTime <- ifelse(is.na(iM0$TTLtime), ttl_default, iM0$TTLtime)

    sN <- parse_wrapped_sweeps(iM0$wrapped_sweeps)

    # extract sweeps
    nwb <- nphys::extractNWB(iMD$nwb_file, sweeps = sN)
    if (is.null(names(nwb$sweeps))) nwb$sweeps <- as.data.frame(nwb$sweeps)
    sweep_names <- names(nwb$sweeps)

    # durations
    dur <- vapply(sweep_names, function(nm) {
      sw <- nwb$sweeps[[nm]]
      sw <- sw[!is.na(sw)]
      sr <- nwb$sampling_rate[[nm]]
      tail(nphys::convertIndex(sw, sr), 1)
    }, numeric(1))

    # APanalysis
    sweep_peaks <- lapply(sweep_names, function(nm) {
      sweep <- nwb$sweeps[[nm]]
      sweep <- sweep[!is.na(sweep)]
      sr <- nwb$sampling_rate[[nm]]
      nphys::APanalysis(sweep, rate = sr)
    })
    names(sweep_peaks) <- sweep_names

    # chunk parse from stimulus_description
    chunk_info <- parse_stim_chunks(iM0$stimulus_description, length(sweep_peaks))
    chunk_names <- chunk_info$parts

    ## Mkke protocol key
    protocol_key <- make_protocol_key(iM0, iMerged, chunk_info, length(sweep_peaks))

    # Build raw trace
    built <- build_raw_trace(sweep_peaks, dur, chunk_names, ttlTime)
    raw_trace <- built$raw_trace
    segment_bounds <- built$segment_bounds

    if (is.null(raw_trace) || nrow(raw_trace) == 0) {
      message("ASP: empty raw trace for ", cell, " | ", protocol_key)
      next
    }

    rec <- recompute_from_raw(raw_trace, iMD, iMerged, iM0)
    tst <- rec$tst
    blMean <- rec$blMean

    # enforce schema in output (guarantee instRate exists)
    tst <- ensure_instRate_percent_change(tst, blMean = blMean)

    # deterministic protocol map
    protocol_map <- list(
      chunk_names = chunk_names,
      baseline_pos = chunk_info$baseline_pos,
      trail_pos = chunk_info$trail_pos,
      ttl_sweep_positions = chunk_info$ttl_pos,
      puff_positions = chunk_info$ttl_pos,
      n_ttl = length(chunk_info$ttl_pos)
    )

    if (isTRUE(plot_it)) {
      graphics::plot(
        tst$time, tst$percent_change,
        main = paste(cell, protocol_key, sep = " - "),
        xlab = "time", ylab = "% Change"
      )
    }

    out_k <- list(
      protocol_key = protocol_key,
      cell = cell,
      md = iMD,
      selected_MD = iMerged,
      pro_row = iM0,
      protocol_map = protocol_map,
      ttlTime = ttlTime,
      sweep_duration = dur,
      sweep_peaks = sweep_peaks,
      segment_bounds = segment_bounds,
      raw_trace = raw_trace,
      spike_puff_output = tst,
      blMean = blMean
    )

    srt$dfs$asp$by_protocol[[protocol_key]] <- out_k

    smry_rows[[protocol_key]] <- data.frame(
      protocol_key = protocol_key,
      stimulus_description = trim1(iM0$stimulus_description),
      n_sweeps = length(sweep_peaks),
      n_ttl = protocol_map$n_ttl,
      ttl_positions = paste(protocol_map$ttl_sweep_positions, collapse = ","),
      wrapped_sweeps = trim1(iM0$wrapped_sweeps),
      TTLtime = ttlTime,
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

  # -----------------------------
  # Header fields + summary table
  # -----------------------------
  srt$dfs$asp$cell <- cell
  srt$dfs$asp$md <- iMD
  srt$dfs$asp$selected_MD <- iMerged
  srt$dfs$asp$protocol_keys <- names(srt$dfs$asp$by_protocol)
  srt$dfs$asp$smry <- if (length(smry_rows)) do.call(rbind, smry_rows) else NULL

  if (isTRUE(save_srt)) {
    saveRDS(srt, file.path(srt$rd, paste0(cell, "-srt.rds")))
  }

  srt$dfs$asp
}
