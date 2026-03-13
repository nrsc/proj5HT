# infer_epochs_from_stimulus <- function(stim, sr,
#                                        on_s = 0.100, off_s = 0.900,
#                                        thr = NULL,
#                                        min_on_s = 0.020) {
#   stopifnot(is.numeric(stim), length(stim) >= 10)
#   stopifnot(is.numeric(sr), length(sr) == 1, sr > 0)
#
#   n <- length(stim)
#   t <- (seq_len(n) - 1) / sr
#
#   # choose threshold if not provided: midpoint between low/high quantiles
#   if (is.null(thr)) {
#     qlo <- stats::quantile(stim, 0.10, na.rm = TRUE, names = FALSE)
#     qhi <- stats::quantile(stim, 0.90, na.rm = TRUE, names = FALSE)
#     thr <- (qlo + qhi) / 2
#   }
#
#   on <- stim > thr
#   on[is.na(on)] <- FALSE
#
#   # run-length encode
#   r <- rle(on)
#   ends <- cumsum(r$lengths)
#   starts <- ends - r$lengths + 1L
#
#   # keep only ON runs above min duration
#   on_runs <- which(r$values)
#   if (length(on_runs) == 0) {
#     # no pulses found
#     return(data.frame(T0 = numeric(0), T1 = numeric(0), subtype = character(0),
#                       epoch = integer(0), pulse = integer(0)))
#   }
#
#   min_on_n <- as.integer(round(min_on_s * sr))
#   keep_runs <- on_runs[r$lengths[on_runs] >= min_on_n]
#   if (length(keep_runs) == 0) keep_runs <- on_runs
#
#   # build segments: Pulse = each ON run; Baseline = from end of ON to next ON start
#   seg <- list()
#   pulse_id <- 0L
#   for (k in seq_along(keep_runs)) {
#     rr <- keep_runs[k]
#     i0 <- starts[rr]
#     i1 <- ends[rr]
#
#     T0p <- t[i0]
#     T1p <- t[i1]
#
#     seg[[length(seg) + 1L]] <- data.frame(
#       T0 = T0p, T1 = T1p, subtype = "Pulse",
#       epoch = 0L, pulse = pulse_id
#     )
#
#     # baseline until next ON
#     if (k < length(keep_runs)) {
#       rr2 <- keep_runs[k + 1L]
#       j0 <- i1 + 1L
#       j1 <- starts[rr2] - 1L
#       if (j1 > j0) {
#         seg[[length(seg) + 1L]] <- data.frame(
#           T0 = t[j0], T1 = t[j1], subtype = "Baseline",
#           epoch = 0L, pulse = pulse_id
#         )
#       }
#     } else {
#       # last baseline until sweep end
#       j0 <- i1 + 1L
#       j1 <- n
#       if (j1 > j0) {
#         seg[[length(seg) + 1L]] <- data.frame(
#           T0 = t[j0], T1 = t[j1], subtype = "Baseline",
#           epoch = 0L, pulse = pulse_id
#         )
#       }
#     }
#
#     pulse_id <- pulse_id + 1L
#   }
#
#   out <- dplyr::bind_rows(seg) %>%
#     dplyr::arrange(T0, T1) %>%
#     dplyr::mutate(
#       # optional: expected duration checks
#       duration_s = T1 - T0
#     )
#
#   out
# }
#
#
# ## Is epoch dataframe missing ##
# is_epochDF_missing <- function(epDF) {
#   # Handles: NULL, 1x1 NA data.frame, or all-NA frames
#   if (is.null(epDF)) return(TRUE)
#   if (!is.data.frame(epDF)) return(TRUE)
#   if (nrow(epDF) == 0) return(TRUE)
#
#   # common bad case: a single NA cell
#   if (nrow(epDF) == 1 && ncol(epDF) == 1) {
#     v <- epDF[[1]][1]
#     if (length(v) == 1 && is.na(v)) return(TRUE)
#   }
#
#   # all entries NA?
#   all_na <- all(is.na(as.matrix(epDF)))
#   if (isTRUE(all_na)) return(TRUE)
#
#   FALSE
# }
#
# get_seg_for_sweep <- function(nwb, nm, sr,
#                               prev_epDF = NULL,
#                               t0_col = "T0", t1_col = "T1",
#                               on_s = 0.100, off_s = 0.900) {
#
#   # Try to fetch epochDF from comments
#   ei <- paste0(nm, ".comment")
#   ep_raw <- NULL
#   if (!is.null(nwb$epochs) && !is.null(nwb$epochs[[ei]]) && "epochDF" %in% names(nwb$epochs[[ei]])) {
#     ep_raw <- nwb$epochs[[ei]][["epochDF"]]
#   }
#
#   epDF <- NULL
#   if (!is.null(ep_raw)) {
#     epDF <- tryCatch(as.data.frame(ep_raw), error = function(e) NULL)
#   }
#
#   # Case A: epochDF is NA -> carry-forward previous valid epochDF
#   if (is_epochDF_missing(epDF)) {
#     if (!is.null(prev_epDF) && !is_epochDF_missing(prev_epDF)) {
#       epDF <- prev_epDF
#     } else {
#       epDF <- NULL
#     }
#   }
#
#   # If we have a usable epDF, build seg via order_epoch_segments
#   if (!is.null(epDF) && !is_epochDF_missing(epDF)) {
#     seg <- order_epoch_segments(epDF, sampling_rate = sr, t0_col = t0_col, t1_col = t1_col)
#     return(list(seg = seg, epDF = epDF, source = "epochDF"))
#   }
#
#   # Case B: no comments / no epochDF -> infer from stimulus trace
#   # ---- You MUST decide where the stimulus vector comes from ----
#   # Common patterns in your NWB extraction:
#   stim <- nwb <- nphys::extractNWB(iMD$nwb_file, stimulus_sweeps = nm, slim = TRUE)
#
#
#   if (is.null(stim)) {
#     stop("No epochDF for sweep ", nm, " and could not find stimulus trace in nwb. ",
#          "Add a stimulus accessor in get_seg_for_sweep().")
#   }
#
#   stim <- as.numeric(stim)
#   seg <- infer_epochs_from_stimulus(stim, sr = sr, on_s = on_s, off_s = off_s)
#
#   # make it compatible with add_pulse_state_time() and (optionally) order_epoch_segments-like output
#   # It already has T0/T1/subtype.
#   seg$i0 <- floor(seg$T0 * sr) + 1L
#   seg$i1 <- floor(seg$T1 * sr)
#   seg$segment_id <- paste0("E", seg$epoch, "_P", seg$pulse, "_", seg$subtype)
#   seg$duration_s <- seg$T1 - seg$T0
#   seg$duration_n <- seg$i1 - seg$i0 + 1L
#
#   list(seg = seg, epDF = NULL, source = "stimulus_infer")
# }
