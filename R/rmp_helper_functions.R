#' RMP helpers (internal)
#'
#' Internal utilities used by RMP building + analysis:
#' parsing, epoch cleaning, pulse-state assignment, stimulus-derived epoching,
#' and raw-trace concatenation.
#'
#' @keywords internal
#' @noRd
NULL

# -----------------------------
# Small utility helpers
# -----------------------------

#' Coalesce NULL to default
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Trim and collapse whitespace
#' @keywords internal
#' @noRd
trim1 <- function(z) trimws(gsub("\\s+", " ", as.character(z)))

#' Safe character conversion: NA -> ""
#' @keywords internal
#' @noRd
safe_chr <- function(z) ifelse(is.na(z), "", as.character(z))

#' Repeat first non-empty element over a vector
#' @keywords internal
#' @noRd
fill_first_nonempty <- function(x, default = NA_character_) {
  x_chr <- trimws(as.character(x))
  x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "nan")] <- NA_character_
  use <- x_chr[!is.na(x_chr)]
  use <- if (length(use) > 0) use[1] else as.character(default)
  rep(use, length(x_chr))
}

#' LOCF with default for leading missings
#' @keywords internal
#' @noRd
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

  if (is.na(x_out[1])) x_out[1] <- default_val
  for (ii in seq_along(x_out)) {
    if (ii == 1) next
    if (is.na(x_out[ii])) x_out[ii] <- x_out[ii - 1]
  }
  x_out[is.na(x_out)] <- default_val
  x_out
}

#' Robust keep column coercion
#' @keywords internal
#' @noRd
as_keep_logical <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))

  x0 <- tolower(trimws(as.character(x)))
  x0[x0 %in% c("", "na", "nan")] <- NA_character_
  out <- rep(FALSE, length(x0))

  out[!is.na(x0) & x0 %in% c("true","t","1","yes","y","keep","k","x")] <- TRUE
  out[!is.na(x0) & x0 %in% c("false","f","0","no","n","drop","d")] <- FALSE

  suppressWarnings({
    xn <- as.numeric(x0)
    out[!is.na(xn)] <- xn[!is.na(xn)] != 0
  })

  out[is.na(x0)] <- FALSE
  out
}

#' Parse wrapped sweeps string into numeric sweep indices (+1 offset)
#' @keywords internal
#' @noRd
parse_wrapped_sweeps <- function(w) {
  w <- as.character(w)[1]
  if (is.na(w) || trimws(w) == "") return(numeric(0))
  w <- gsub("\\[|\\]", "", w)
  sp <- trimws(unlist(strsplit(w, ",")))
  sp <- sp[sp != ""]
  suppressWarnings(as.numeric(sp)) + 1
}

# -----------------------------
# RMP-specific helpers
# -----------------------------

#' Parse stim chunks for RMP protocols
#'
#' Splits `stimulus_description` by `+` and aligns to `n_sweeps`.
#' Baseline chunk is expected as `tpRMPbl`, puff chunk `tpRMP_puff`, and `trail`.
#'
#' @keywords internal
#' @noRd
parse_stim_chunks_rmp <- function(stim_desc, n_sweeps) {
  parts0 <- unlist(strsplit(as.character(stim_desc), "\\+"))
  parts0 <- trimws(parts0)
  parts0 <- parts0[parts0 != ""]

  if (length(parts0) != n_sweeps) {
    if (length(parts0) == 0) parts0 <- rep("unknown", n_sweeps)
    if (length(parts0) < n_sweeps) parts0 <- c(parts0, rep("trail", n_sweeps - length(parts0)))
    if (length(parts0) > n_sweeps) parts0 <- parts0[seq_len(n_sweeps)]
  }

  list(
    parts = parts0,
    puff_pos = which(parts0 == "tpRMP_puff"),
    baseline_pos = which(parts0 == "tpRMPbl"),
    trail_pos = which(parts0 == "trail")
  )
}

#' Sweep vector -> tpRMP dataframe with `epoch_t`
#' @keywords internal
#' @noRd
sweep_to_tpRMP_df <- function(sweep_vec, sr) {
  sweep_vec <- sweep_vec[!is.na(sweep_vec)]
  n <- length(sweep_vec)
  if (n == 0) return(data.frame(i = integer(0), epoch_t = numeric(0), tpRMP = numeric(0)))
  i <- seq_len(n)
  t <- (i - 1) / sr
  data.frame(i = i, epoch_t = t, tpRMP = as.numeric(sweep_vec))
}

#' Clean an epochDF (drop header rows, coerce numeric times)
#'
#' Expects columns like T0/T1 plus V3/V4/V5 (from your comments epochDF format).
#'
#' @keywords internal
#' @noRd
clean_epochDF <- function(epDF, t0_col = "T0", t1_col = "T1") {
  if (is.null(epDF) || !is.data.frame(epDF) || nrow(epDF) == 0) return(NULL)

  need <- c(t0_col, t1_col, "V3", "V4", "V5")
  if (!all(need %in% names(epDF))) return(NULL)

  to_num <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

  ep <- epDF
  ep[[t0_col]] <- to_num(ep[[t0_col]])
  ep[[t1_col]] <- to_num(ep[[t1_col]])

  ep <- ep[is.finite(ep[[t0_col]]) & is.finite(ep[[t1_col]]) & ep[[t1_col]] > ep[[t0_col]], , drop = FALSE]
  if (nrow(ep) == 0) return(NULL)

  ep
}

#' Add pulse_state using time column and segment (T0/T1) bounds
#'
#' Marks "pulse_on" where seg subtype is Pulse and `epoch_t` is within [T0, T1].
#'
#' @keywords internal
#' @noRd
add_pulse_state_time <- function(df, seg,
                                 df_time_col = "epoch_t",
                                 seg_t0_col = "T0",
                                 seg_t1_col = "T1",
                                 seg_subtype_col = "subtype",
                                 state_col = "pulse_state") {
  stopifnot(is.data.frame(df), is.data.frame(seg))
  stopifnot(df_time_col %in% names(df))
  stopifnot(all(c(seg_t0_col, seg_t1_col, seg_subtype_col) %in% names(seg)))

  t <- as.numeric(df[[df_time_col]])
  state <- rep("pulse_off", length(t))

  seg2 <- seg
  seg2[[seg_subtype_col]] <- trimws(as.character(seg2[[seg_subtype_col]]))

  seg_on <- seg2[seg2[[seg_subtype_col]] == "Pulse" &
                   is.finite(seg2[[seg_t0_col]]) & is.finite(seg2[[seg_t1_col]]), , drop = FALSE]

  if (nrow(seg_on) > 0) {
    for (k in seq_len(nrow(seg_on))) {
      hit <- t >= seg_on[[seg_t0_col]][k] & t <= seg_on[[seg_t1_col]][k]
      state[hit] <- "pulse_on"
    }
  }

  df[[state_col]] <- state
  df
}

#' Make epochDF from a 2-level stimulus vector by detecting pulse ON segments
#'
#' @keywords internal
#' @noRd
make_epochDF_from_stimulus <- function(stim, sr,
                                       on_label = "Pulse",
                                       off_label = "Baseline",
                                       thr = NULL,
                                       min_on_s = 0.02,
                                       expected_on_s = 0.100,
                                       pulse_polarity = c("auto","high","low")) {

  pulse_polarity <- match.arg(pulse_polarity)
  stopifnot(is.numeric(stim), length(stim) >= 10)
  stopifnot(is.numeric(sr), length(sr) == 1, sr > 0)

  n <- length(stim)
  t <- (seq_len(n) - 1) / sr

  if (is.null(thr)) {
    qlo <- stats::quantile(stim, 0.10, na.rm = TRUE, names = FALSE)
    qhi <- stats::quantile(stim, 0.90, na.rm = TRUE, names = FALSE)
    thr <- (qlo + qhi) / 2
  }

  # candidate ON masks
  on_high <- stim > thr     # "high" state is ON (e.g. 0 vs -100)
  on_low  <- stim < thr     # "low"  state is ON (e.g. -100 vs 0)
  on_high[is.na(on_high)] <- FALSE
  on_low[is.na(on_low)]   <- FALSE

  choose_mask <- function(on_mask) {
    r <- rle(on_mask)
    if (!any(r$values)) return(list(score = Inf, mask = on_mask))
    on_lens <- r$lengths[r$values]
    # score: how close is the median ON run to expected_on_s?
    med_s <- stats::median(on_lens) / sr
    list(score = abs(med_s - expected_on_s), mask = on_mask)
  }

  if (pulse_polarity == "high") {
    on <- on_high
  } else if (pulse_polarity == "low") {
    on <- on_low
  } else {
    a <- choose_mask(on_high)
    b <- choose_mask(on_low)
    on <- if (a$score <= b$score) a$mask else b$mask
  }

  r <- rle(on)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1L

  on_runs <- which(r$values)
  min_on_n <- as.integer(round(min_on_s * sr))
  on_runs <- on_runs[r$lengths[on_runs] >= min_on_n]
  if (length(on_runs) == 0) return(NULL)

  rows <- list()
  pulse_id <- 0L

  for (k in seq_along(on_runs)) {
    rr <- on_runs[k]
    i0 <- starts[rr]
    i1 <- ends[rr]

    T0p <- t[i0]
    T1p <- t[i1]

    rows[[length(rows) + 1L]] <- data.frame(
      T0 = sprintf("%.8f", T0p),
      T1 = sprintf("%.8f", T1p),
      Type = "Type=Epoch",
      Subtype = "EpochType=Pulse Train",
      V1 = "Epoch=0",
      V2 = NA_character_,
      V3 = paste0("Pulse=", pulse_id),
      V4 = paste0("SubType=", on_label),
      V5 = paste0("ShortName=E0_PT_P", pulse_id, "_P"),
      stringsAsFactors = FALSE
    )

    # baseline after pulse until next pulse or end
    if (k < length(on_runs)) {
      rr2 <- on_runs[k + 1L]
      j0 <- i1 + 1L
      j1 <- starts[rr2] - 1L
    } else {
      j0 <- i1 + 1L
      j1 <- n
    }

    if (j1 > j0) {
      rows[[length(rows) + 1L]] <- data.frame(
        T0 = sprintf("%.8f", t[j0]),
        T1 = sprintf("%.8f", t[j1]),
        Type = "Type=Epoch",
        Subtype = "EpochType=Pulse Train",
        V1 = "Epoch=0",
        V2 = NA_character_,
        V3 = paste0("Pulse=", pulse_id),
        V4 = paste0("SubType=", off_label),
        V5 = paste0("ShortName=E0_PT_P", pulse_id, "_B"),
        stringsAsFactors = FALSE
      )
    }

    pulse_id <- pulse_id + 1L
  }

  epDF <- dplyr::bind_rows(rows)

  header <- data.frame(
    T0 = c("HS#0", "Epochs"),
    T1 = c(NA, NA),
    Type = c(NA, NA),
    Subtype = c(NA, NA),
    V1 = c(NA, NA),
    V2 = c(NA, NA),
    V3 = c(NA, NA),
    V4 = c(NA, NA),
    V5 = c(NA, NA),
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(header, epDF)
}



#' Build concatenated raw trace for RMP protocol set
#'
#' Concatenates sweeps with a global time axis aligned so the first tpRMP_puff
#' sweep is near zero, then shifts by `puffTime` so puff onset lands at ~0.
#'
#' @keywords internal
#' @noRd
build_raw_trace_rmp <- function(sweep_objs, dur, chunk_names, puffTime) {
  sweep_dfs <- lapply(sweep_objs, `[[`, "df")

  cum_shift <- c(0, cumsum(as.numeric(dur)))[seq_along(sweep_dfs)]

  first_puff_idx <- which(chunk_names == "tpRMP_puff")[1]
  if (is.na(first_puff_idx)) first_puff_idx <- 1L
  zero_offset <- cum_shift[first_puff_idx]

  segment_bounds <- data.frame(
    sweep = seq_along(sweep_dfs),
    chunk = chunk_names,
    start_time = NA_real_,
    end_time = NA_real_,
    puff_onset_time = NA_real_,
    stringsAsFactors = FALSE
  )

  out_list <- vector("list", length(sweep_dfs))

  for (i in seq_along(sweep_dfs)) {
    d <- sweep_dfs[[i]]
    t_pre <- d$epoch_t + cum_shift[i] - zero_offset
    t <- t_pre - puffTime

    out_list[[i]] <- data.frame(
      time = t,
      tpRMP = d$tpRMP,
      pulse_state = d$pulse_state %||% NA_character_,
      chunk = chunk_names[i],
      sweep = i,
      stringsAsFactors = FALSE
    )

    segment_bounds$start_time[i] <- min(t, na.rm = TRUE)
    segment_bounds$end_time[i]   <- max(t, na.rm = TRUE)

    if (chunk_names[i] == "tpRMP_puff" && nrow(d) > 0) {
      segment_bounds$puff_onset_time[i] <- 0
    }
  }

  raw <- do.call(rbind, out_list)
  list(raw_trace = raw, segment_bounds = segment_bounds)
}

#' Recompute aligned analysis table and baseline mean from raw_trace
#'
#' Applies an outlier filter on tpRMP and computes `protocol` and `blMean`.
#'
#' @keywords internal
#' @noRd
recompute_from_raw_rmp <- function(raw_trace, outlier_sd = 5, bl_window = c(-5, 0)) {
  raw_trace <- raw_trace[!is.na(raw_trace$time) & !is.na(raw_trace$tpRMP), , drop = FALSE]

  tst <- raw_trace
  mu  <- mean(tst$tpRMP, na.rm = TRUE)
  sig <- stats::sd(tst$tpRMP, na.rm = TRUE)
  if (is.finite(sig) && sig > 0) {
    keep <- abs(tst$tpRMP - mu) <= outlier_sd * sig
    if (any(!keep, na.rm = TRUE)) tst <- tst[keep, , drop = FALSE]
  }

  tst$protocol <- ifelse(tst$time <= 0, "baseline", "stimulus")

  bl_keep <- tst$time > bl_window[1] & tst$time <= bl_window[2]
  blMean  <- mean(tst$tpRMP[bl_keep], na.rm = TRUE)

  list(tst = tst, blMean = blMean)
}

ensure_asp <- function(x,
                       mMD = NULL,
                       mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                       mMD_sheet = "Master_Macaque",
                       ttl_default = 5,
                       outlier_sd = 5,
                       bl_window = c(-5, 0),
                       save_srt = FALSE) {

  # accept cell_name, srt, or asp
  if (is.character(x)) {
    srt <- load5HT(x, tag = "-srt.rds")
  } else if (is.list(x) && !is.null(x$by_protocol) && !is.null(x$protocol_keys)) {
    # looks like an asp object already
    return(list(srt = NULL, asp = x))
  } else {
    srt <- x
  }

  # already built?
  if (!is.null(srt$dfs$asp) && is.list(srt$dfs$asp) && length(srt$dfs$asp$by_protocol) > 0) {
    return(list(srt = srt, asp = srt$dfs$asp))
  }

  asp <- build_asp(
    srt,
    mMD = mMD,
    mMD_path = mMD_path,
    mMD_sheet = mMD_sheet,
    ttl_default = ttl_default,
    outlier_sd = outlier_sd,
    bl_window = bl_window,
    save_srt = save_srt
  )

  # build_asp returns asp list; also stored in srt if save_srt TRUE
  srt$dfs$asp <- asp
  list(srt = srt, asp = asp)
}


pick_asp_protocols <- function(asp,
                               rank = NULL,
                               exp = NULL,
                               stim_grepl = NULL,
                               keep_only = TRUE) {

  keys <- names(asp$by_protocol)

  meta <- dplyr::bind_rows(lapply(keys, function(k) {
    pr <- asp$by_protocol[[k]]$pro_row
    data.frame(
      protocol_key = k,
      Rank = if ("Rank" %in% names(pr)) as.character(pr$Rank) else NA_character_,
      Exp  = if ("Exp"  %in% names(pr)) as.character(pr$Exp)  else NA_character_,
      stimulus_description = if ("stimulus_description" %in% names(pr)) as.character(pr$stimulus_description) else NA_character_,
      keep = if ("keep" %in% names(pr)) isTRUE(pr$keep) else NA,
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(rank)) {
    meta <- meta[dplyr::if_else(is.na(meta$Rank), FALSE, meta$Rank %in% rank), , drop = FALSE]
  }
  if (!is.null(exp)) {
    meta <- meta[dplyr::if_else(is.na(meta$Exp), FALSE, meta$Exp %in% exp), , drop = FALSE]
  }
  if (!is.null(stim_grepl)) {
    meta <- meta[grepl(stim_grepl, meta$stimulus_description), , drop = FALSE]
  }
  if (isTRUE(keep_only) && "keep" %in% names(meta)) {
    meta <- meta[isTRUE(meta$keep) | is.na(meta$keep), , drop = FALSE]
  }

  meta
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

stop_if_not_srt <- function(srt) {
  if (!is.list(srt)) stop("Expected srt to be a list/nnest.")
  if (is.null(srt$cell)) stop("srt$cell is missing.")
  if (is.null(srt$rd))   stop("srt$rd is missing (needed for saving).")
  invisible(TRUE)
}

ensure_asp_in_srt <- function(srt,
                              mMD = NULL,
                              mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                              mMD_sheet = "Master_Macaque",
                              ttl_default = 5,
                              outlier_sd = 5,
                              bl_window = c(-5, 0),
                              save_srt = TRUE) {
  stop_if_not_srt(srt)

  has_asp <- !is.null(srt$dfs$asp) &&
    is.list(srt$dfs$asp) &&
    !is.null(srt$dfs$asp$by_protocol) &&
    length(srt$dfs$asp$by_protocol) > 0

  if (isTRUE(has_asp)) return(srt)

  asp <- build_asp(
    srt,
    mMD = mMD,
    mMD_path = mMD_path,
    mMD_sheet = mMD_sheet,
    ttl_default = ttl_default,
    outlier_sd = outlier_sd,
    bl_window = bl_window,
    save_srt = save_srt
  )

  # build_asp already stores at srt$dfs$asp if save_srt TRUE (in your code it returns asp only)
  # so we defensively set it here too:
  srt$dfs$asp <- asp

  if (isTRUE(save_srt)) {
    saveRDS(srt, file.path(srt$rd, paste0(srt$cell, "-srt.rds")))
  }

  srt
}

order_epoch_segments <- function(epochDF, sampling_rate,
                                 keep_subtypes = c("Baseline", "Pulse"),
                                 t0_col = NULL, t1_col = NULL) {

  stopifnot(is.data.frame(epochDF))
  stopifnot(is.numeric(sampling_rate), length(sampling_rate) == 1, sampling_rate > 0)

  nm <- names(epochDF)

  # ---- guess time columns if not provided ----
  if (is.null(t0_col) || is.null(t1_col)) {
    candidates <- list(
      c("T0","T1"),
      c("t0","t1"),
      c("start_time","stop_time"),
      c("start_time_s","stop_time_s"),
      c("start","end"),
      c("time_start","time_end"),
      c("t_start","t_end")
    )

    found <- FALSE
    for (pair in candidates) {
      if (all(pair %in% nm)) {
        t0_col <- pair[1]
        t1_col <- pair[2]
        found <- TRUE
        break
      }
    }

    if (!found) {
      stop(
        "Couldn't find time columns. Pass t0_col/t1_col explicitly. ",
        "Available columns: ", paste(nm, collapse = ", ")
      )
    }
  }

  # Helper: pull "Key=Value" tokens out of free text
  parse_kv <- function(x, key) {
    out <- sub(paste0(".*\\b", key, "=([^ ]+).*"), "\\1", x)
    out[out == x] <- NA_character_
    out
  }

  df <- epochDF

  # Ensure numeric times
  df[[t0_col]] <- as.numeric(df[[t0_col]])
  df[[t1_col]] <- as.numeric(df[[t1_col]])

  # Parse Epoch / Pulse / SubType from V3/V4/V5 if present
  txt <- paste(df$V3 %||% "", df$V4 %||% "", df$V5 %||% "")

  df$epoch   <- suppressWarnings(as.integer(parse_kv(txt, "Epoch")))
  df$pulse   <- suppressWarnings(as.integer(parse_kv(txt, "Pulse")))
  df$subtype <- parse_kv(txt, "SubType")

  # Keep only the rows that actually define segments
  df <- df |>
    dplyr::filter(!is.na(.data[[t0_col]]), !is.na(.data[[t1_col]])) |>
    dplyr::filter(.data[[t1_col]] > .data[[t0_col]]) |>
    dplyr::filter(!is.na(subtype), subtype %in% keep_subtypes) |>
    dplyr::arrange(.data[[t0_col]], .data[[t1_col]], epoch, pulse, subtype) |>
    dplyr::mutate(
      i0 = floor(.data[[t0_col]] * sampling_rate) + 1L,
      i1 = floor(.data[[t1_col]] * sampling_rate),
      i0 = pmax(i0, 1L),
      i1 = pmax(i1, i0),
      segment_id = paste0("E", epoch, "_P", pulse, "_", subtype),
      duration_s = .data[[t1_col]] - .data[[t0_col]],
      duration_n = i1 - i0 + 1L
    )

  df
}

get_seg_for_sweep <- function(nwb, nm, sr,
                              prev_epDF = NULL,
                              t0_col = "T0", t1_col = "T1",
                              on_s = 0.100, off_s = 0.900) {

  # Try to fetch epochDF from comments
  ei <- paste0(nm, ".comment")
  ep_raw <- NULL
  if (!is.null(nwb$epochs) && !is.null(nwb$epochs[[ei]]) && "epochDF" %in% names(nwb$epochs[[ei]])) {
    ep_raw <- nwb$epochs[[ei]][["epochDF"]]
  }

  epDF <- NULL
  if (!is.null(ep_raw)) {
    epDF <- tryCatch(as.data.frame(ep_raw), error = function(e) NULL)
  }

  # Case A: epochDF is NA -> carry-forward previous valid epochDF
  if (is_epochDF_missing(epDF)) {
    if (!is.null(prev_epDF) && !is_epochDF_missing(prev_epDF)) {
      epDF <- prev_epDF
    } else {
      epDF <- NULL
    }
  }

  # If we have a usable epDF, build seg via order_epoch_segments
  if (!is.null(epDF) && !is_epochDF_missing(epDF)) {
    seg <- order_epoch_segments(epDF, sampling_rate = sr, t0_col = t0_col, t1_col = t1_col)
    return(list(seg = seg, epDF = epDF, source = "epochDF"))
  }

  # Case B: no comments / no epochDF -> infer from stimulus trace
  # ---- You MUST decide where the stimulus vector comes from ----
  # Common patterns in your NWB extraction:
  stim <- nwb <- nphys::extractNWB(iMD$nwb_file, stimulus_sweeps = nm, slim = TRUE)


  if (is.null(stim)) {
    stop("No epochDF for sweep ", nm, " and could not find stimulus trace in nwb. ",
         "Add a stimulus accessor in get_seg_for_sweep().")
  }

  stim <- as.numeric(stim)
  seg <- infer_epochs_from_stimulus(stim, sr = sr, on_s = on_s, off_s = off_s)

  # make it compatible with add_pulse_state_time() and (optionally) order_epoch_segments-like output
  # It already has T0/T1/subtype.
  seg$i0 <- floor(seg$T0 * sr) + 1L
  seg$i1 <- floor(seg$T1 * sr)
  seg$segment_id <- paste0("E", seg$epoch, "_P", seg$pulse, "_", seg$subtype)
  seg$duration_s <- seg$T1 - seg$T0
  seg$duration_n <- seg$i1 - seg$i0 + 1L

  list(seg = seg, epDF = NULL, source = "stimulus_infer")
}




#' Construct a deterministic protocol key for RMP protocol sets
#'
#' Builds a human-readable, stable key encoding:
#' protocol type, subclass, sweep count, puff count, sweep list, puff time,
#' bath drug, puff drug, and experiment label.
#'
#' @param iM0 One-row metadata (Master sheet row) for this protocol.
#' @param chunk_info List from `parse_stim_chunks_rmp()`.
#' @param n_sweeps Number of sweeps in this protocol set.
#' @param puff_default Default puff time if TTLtime missing/NA.
#'
#' @return Character scalar protocol key (sanitized).
#' @keywords internal
#' @noRd
make_protocol_key_rmp <- function(iM0, iMD, chunk_info, n_sweeps, puff_default = 5) {

  stim_clean <- trim1(iM0$stimulus_description[1])
  typ <- if (grepl("tpRMP_puff", stim_clean, fixed = TRUE)) "tpRMP" else "other"

  puff_t <- iM0$TTLtime[1]
  if (!is.finite(suppressWarnings(as.numeric(puff_t))) || is.na(puff_t)) puff_t <- puff_default
  puff_t <- as.numeric(puff_t)

  key <- paste(
    typ,
    paste0("subc=", iM0$assigned_subclass[1]),
    paste0("nSweep=", n_sweeps),
    paste0("nPuff=", length(chunk_info$puff_pos)),
    paste0("sweeps=", trim1(iM0$wrapped_sweeps[1])),
    paste0("puffTime=", puff_t),
    paste0("bath=", safe_chr(iMD$bath_Drug[1])),
    paste0("puffdrug=", safe_chr(iMD$puff_drug[1])),
    paste0("exp=", safe_chr(iM0$Exp[1])),
    sep = "|"
  )

  # sanitize: keep alnum + a few separators, convert the rest to underscores
  gsub("[^A-Za-z0-9|=_\\-]+", "_", key)
}

is_epochDF_usable <- function(epDF, t0_col = "T0", t1_col = "T1") {
  ep <- clean_epochDF(epDF, t0_col = t0_col, t1_col = t1_col)
  !is.null(ep) && nrow(ep) > 0
}


get_stim_for_sweep <- function(nwb, nm) {
  stim <- nphys::extractNWB(nwb$nwb_file, stimulus_sweeps = nm, slim = TRUE) # Try common places first:
  if (!is.null(stim)) return(as.numeric(stim))
  NULL
}
#
# m_idx <- which(as.character(mMD$cell_name) == cell)
# if (length(m_idx) == 0) {
#   message("No Master_UI_data match for cell: ", cell)
#   return(NULL)
# }
# if (length(m_idx) > 1) {
#   #message("Multiple Master_UI_data matches for cell ", cell, " (keeping all).")
# }
# iMerged <- mMD[m_idx, , drop = FALSE]
# iMerged$assigned_subclass = iMerged$assigned_subclass[1]
#
# # columns you want to compare
# cols_keep <- c(1:3, 5, 6, 7,13)
#
# same_as_asp <- !is.null(srt$dfs$asp) &&
#   identical(iMerged[, cols_keep], srt$dfs$asp$selected_MD[, cols_keep])
#
# same_as_rmp <- !is.null(srt$dfs$rmp) &&
#   identical(iMerged[, cols_keep], srt$dfs$rmp$selected_MD[, cols_keep])
#
# # skip only if ALL available comparisons are identical
# if (same_as_asp || same_as_rmp) {
#   # If both exist, require both to be identical to skip
#   if ((!is.null(srt$dfs$asp) && !is.null(srt$dfs$rmp) && same_as_asp && same_as_rmp) ||
#       (!is.null(srt$dfs$asp) &&  is.null(srt$dfs$rmp) && same_as_asp) ||
#       ( is.null(srt$dfs$asp) && !is.null(srt$dfs$rmp) && same_as_rmp)) {
#     message("no changes to iMerged")
#     next
#   }
# }

