# # build_rmp_helper_functions
#
# # -----------------------------
# # Helper functions
# # -----------------------------
# trim1 <- function(z) trimws(gsub("\\s+", " ", as.character(z)))
# safe_chr <- function(z) ifelse(is.na(z), "", as.character(z))
#
# # Repeat first non-empty value for the whole column (per-cell constants)
# fill_first_nonempty <- function(x, default = NA_character_) {
#   x_chr <- trimws(as.character(x))
#   x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "nan")] <- NA_character_
#   use <- x_chr[!is.na(x_chr)]
#   use <- if (length(use) > 0) use[1] else as.character(default)
#   rep(use, length(x_chr))
# }
#
# # LOCF with default for leading missings
# fill_locf <- function(x, default, type = c("character", "numeric")) {
#   type <- match.arg(type)
#
#   x_chr <- trimws(as.character(x))
#   x_chr[x_chr == "" | tolower(x_chr) %in% c("na", "nan")] <- NA_character_
#
#   if (type == "numeric") {
#     x_out <- suppressWarnings(as.numeric(x_chr))
#     default_val <- as.numeric(default)
#   } else {
#     x_out <- x_chr
#     default_val <- as.character(default)
#   }
#
#   if (length(x_out) == 0L) return(x_out)
#
#   if (is.na(x_out[1])) x_out[1] <- default_val
#   for (ii in seq_along(x_out)) {
#     if (ii == 1) next
#     if (is.na(x_out[ii])) x_out[ii] <- x_out[ii - 1]
#   }
#   x_out[is.na(x_out)] <- default_val
#   x_out
# }
#
# as_keep_logical <- function(x) {
#   if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
#
#   x0 <- tolower(trimws(as.character(x)))
#   x0[x0 %in% c("", "na", "nan")] <- NA_character_
#   out <- rep(FALSE, length(x0))
#
#   out[!is.na(x0) & x0 %in% c("true","t","1","yes","y","keep","k","x")] <- TRUE
#   out[!is.na(x0) & x0 %in% c("false","f","0","no","n","drop","d")] <- FALSE
#
#   suppressWarnings({
#     xn <- as.numeric(x0)
#     out[!is.na(xn)] <- xn[!is.na(xn)] != 0
#   })
#
#   out[is.na(x0)] <- FALSE
#   out
# }
#
# parse_wrapped_sweeps <- function(w) {
#   w <- as.character(w)[1]
#   if (!is.finite(suppressWarnings(as.numeric(w))) && (is.na(w) || trimws(w) == "")) return(numeric(0))
#   w <- gsub("\\[|\\]", "", w)
#   sp <- trimws(unlist(strsplit(w, ",")))
#   sp <- sp[sp != ""]
#   suppressWarnings(as.numeric(sp)) + 1
# }
#
# # Define chunk names + puff sweep positions from stim_desc
# # Baseline chunk is `tpRMP_bl`, puff chunk is `tpRMP_puff`, and keep `trail`.
# parse_stim_chunks_rmp <- function(stim_desc, n_sweeps) {
#   parts0 <- unlist(strsplit(as.character(stim_desc), "\\+"))
#   parts0 <- trimws(parts0)
#   parts0 <- parts0[parts0 != ""]
#
#   if (length(parts0) != n_sweeps) {
#     if (length(parts0) == 0) parts0 <- rep("unknown", n_sweeps)
#     if (length(parts0) < n_sweeps) parts0 <- c(parts0, rep("trail", n_sweeps - length(parts0)))
#     if (length(parts0) > n_sweeps) parts0 <- parts0[seq_len(n_sweeps)]
#   }
#
#   puff_pos <- which(parts0 == "tpRMP_puff")
#   list(
#     parts = parts0,
#     puff_pos = puff_pos,
#     baseline_pos = which(parts0 == "tpRMP_bl"),
#     trail_pos = which(parts0 == "trail")
#   )
# }
#
# # Turn a sweep vector into a time series DF (tpRMP == voltage trace)
# sweep_to_tpRMP_df <- function(sweep_vec, sr) {
#   sweep_vec <- sweep_vec[!is.na(sweep_vec)]
#   n <- length(sweep_vec)
#   if (n == 0) return(data.frame(i = integer(0), epoch_t = numeric(0), tpRMP = numeric(0)))
#   i <- seq_len(n)
#   t <- (i - 1) / sr
#   data.frame(i = i, epoch_t = t, tpRMP = as.numeric(sweep_vec))
# }
#
#
# order_epoch_segments <- function(epochDF, sampling_rate,
#                                  keep_subtypes = c("Baseline", "Pulse"),
#                                  t0_col = NULL, t1_col = NULL) {
#
#   stopifnot(is.data.frame(epochDF))
#   stopifnot(is.numeric(sampling_rate), length(sampling_rate) == 1, sampling_rate > 0)
#
#   nm <- names(epochDF)
#
#   # ---- guess time columns if not provided ----
#   if (is.null(t0_col) || is.null(t1_col)) {
#     candidates <- list(
#       c("T0","T1"),
#       c("t0","t1"),
#       c("start_time","stop_time"),
#       c("start_time_s","stop_time_s"),
#       c("start","end"),
#       c("time_start","time_end"),
#       c("t_start","t_end")
#     )
#
#     found <- FALSE
#     for (pair in candidates) {
#       if (all(pair %in% nm)) {
#         t0_col <- pair[1]
#         t1_col <- pair[2]
#         found <- TRUE
#         break
#       }
#     }
#
#     if (!found) {
#       stop(
#         "Couldn't find time columns. Pass t0_col/t1_col explicitly. ",
#         "Available columns: ", paste(nm, collapse = ", ")
#       )
#     }
#   }
#
#   # Helper: pull "Key=Value" tokens out of free text
#   parse_kv <- function(x, key) {
#     out <- sub(paste0(".*\\b", key, "=([^ ]+).*"), "\\1", x)
#     out[out == x] <- NA_character_
#     out
#   }
#
#   df <- epochDF
#
#   # Ensure numeric times
#   df[[t0_col]] <- as.numeric(df[[t0_col]])
#   df[[t1_col]] <- as.numeric(df[[t1_col]])
#
#   # Parse Epoch / Pulse / SubType from V3/V4/V5 if present
#   txt <- paste(df$V3 %||% "", df$V4 %||% "", df$V5 %||% "")
#
#   df$epoch   <- suppressWarnings(as.integer(parse_kv(txt, "Epoch")))
#   df$pulse   <- suppressWarnings(as.integer(parse_kv(txt, "Pulse")))
#   df$subtype <- parse_kv(txt, "SubType")
#
#   # Keep only the rows that actually define segments
#   df <- df |>
#     dplyr::filter(!is.na(.data[[t0_col]]), !is.na(.data[[t1_col]])) |>
#     dplyr::filter(.data[[t1_col]] > .data[[t0_col]]) |>
#     dplyr::filter(!is.na(subtype), subtype %in% keep_subtypes) |>
#     dplyr::arrange(.data[[t0_col]], .data[[t1_col]], epoch, pulse, subtype) |>
#     dplyr::mutate(
#       i0 = floor(.data[[t0_col]] * sampling_rate) + 1L,
#       i1 = floor(.data[[t1_col]] * sampling_rate),
#       i0 = pmax(i0, 1L),
#       i1 = pmax(i1, i0),
#       segment_id = paste0("E", epoch, "_P", pulse, "_", subtype),
#       duration_s = .data[[t1_col]] - .data[[t0_col]],
#       duration_n = i1 - i0 + 1L
#     )
#
#   df
# }
#
# add_pulse_state <- function(df, seg,
#                             index_col = NULL,
#                             time_col = NULL,
#                             sampling_rate = NULL,
#                             state_col = "pulse_state",
#                             id_col = NULL) {
#
#   stopifnot(is.data.frame(df), is.data.frame(seg))
#   stopifnot(all(c("i0","i1","subtype") %in% names(seg)))
#
#   # Decide how we're indexing df: by i (preferred) or by time
#   if (is.null(index_col)) {
#     idx_candidates <- intersect(c("i","idx","sample","sample_idx"), names(df))
#     if (length(idx_candidates) > 0) index_col <- idx_candidates[1]
#   }
#   if (is.null(time_col)) {
#     time_candidates <- intersect(c("time","t","seconds"), names(df))
#     if (length(time_candidates) > 0) time_col <- time_candidates[1]
#   }
#
#   # Build index vector
#   if (!is.null(index_col)) {
#     idx <- as.integer(df[[index_col]])
#   } else if (!is.null(time_col)) {
#     stopifnot(is.numeric(sampling_rate), length(sampling_rate) == 1, sampling_rate > 0)
#     idx <- floor(as.numeric(df[[time_col]]) * sampling_rate) + 1L
#   } else {
#     # FINAL fallback: assume each row is one sample in order
#     idx <- seq_len(nrow(df))
#   }
#
#   n <- length(idx)
#   state <- rep("pulse_off", n)
#
#   # only Pulse segments define "on"
#   seg_on <- seg[seg$subtype == "Pulse" & !is.na(seg$i0) & !is.na(seg$i1), , drop = FALSE]
#   if (nrow(seg_on) > 0) {
#     for (k in seq_len(nrow(seg_on))) {
#       i0 <- seg_on$i0[k]
#       i1 <- seg_on$i1[k]
#       hit <- idx >= i0 & idx <= i1
#       state[hit] <- "pulse_on"
#     }
#   }
#
#   df[[state_col]] <- state
#
#   if (!is.null(id_col)) {
#     df[[id_col]] <- NA_character_
#     if (nrow(seg_on) > 0) {
#       for (k in seq_len(nrow(seg_on))) {
#         i0 <- seg_on$i0[k]
#         i1 <- seg_on$i1[k]
#         hit <- idx >= i0 & idx <= i1
#         df[[id_col]][hit] <- paste0("E", seg_on$epoch[k], "_P", seg_on$pulse[k])
#       }
#     }
#   }
#
#   df
# }
#
# recompute_from_raw_rmp <- function(raw_trace, iMD, iMerged, iM0) {
#   raw_trace <- raw_trace[!is.na(raw_trace$time) & !is.na(raw_trace$tpRMP), , drop = FALSE]
#
#   # outlier filter on tpRMP
#   tst <- raw_trace
#   mu  <- mean(tst$tpRMP, na.rm = TRUE)
#   sig <- stats::sd(tst$tpRMP, na.rm = TRUE)
#   if (is.finite(sig) && sig > 0) {
#     keep <- abs(tst$tpRMP - mu) <= outlier_sd * sig
#     if (any(!keep, na.rm = TRUE)) tst <- tst[keep, , drop = FALSE]
#   }
#
#   # protocol flag (baseline vs stimulus based on aligned time)
#   tst$protocol <- ifelse(tst$time <= 0, "baseline", "stimulus")
#
#   # baseline mean from aligned time window
#   bl_keep <- tst$time > bl_window[1] & tst$time <= bl_window[2]
#   blMean  <- mean(tst$tpRMP[bl_keep], na.rm = TRUE)
#
#   if (!is.finite(blMean)) {
#     blMean <- NA_real_
#     tst$delta_mV <- NA_real_
#     tst$percent_of_baseline <- NA_real_
#   } else {
#     tst$delta_mV <- tst$tpRMP - blMean
#     if (abs(blMean) < 1e-6) {
#       tst$percent_of_baseline <- NA_real_
#     } else {
#       tst$percent_of_baseline <- (tst$tpRMP / blMean) * 100
#     }
#   }
#
#   # assigned_subclass: per-cell constant from iMerged
#   assigned_sub <- iMerged$assigned_subclass
#   assigned_sub <- assigned_sub[!is.na(assigned_sub)]
#   assigned_sub <- if (length(assigned_sub) == 0) NA_character_ else as.character(assigned_sub[1])
#
#   # metadata columns (mirrors your ASP output)
#   tst$cell_name <- iMD$cell_name
#   tst$predicted_subclass <- iMD$predicted_subclass
#   tst$Species <- iMD$Species
#   tst$blockers <- iMD$Fast_blockers
#   tst$bath <- iM0$Bath_Drug
#   tst$puff <- iM0$Puff_drug
#   tst$DFP <- iMD$Depth_from_pia
#   tst$LIMs_depth <- iMD$LIMS_depth
#   tst$assigned_subclass <- rep(assigned_sub, nrow(tst))
#   tst$expCon <- safe_chr(as.character(iM0$Exp[1]))
#   tst$response <- iM0$Effect
#   tst$keep <- as_keep_logical(iM0$keep)[1]
#
#   list(tst = tst, blMean = blMean)
# }
#
# make_protocol_key_rmp <- function(iM0, chunk_info, n_sweeps) {
#   stim_clean <- trim1(iM0$stimulus_description[1])
#   typ <- if (grepl("tpRMP_puff", stim_clean)) "tpRMP" else "other"
#   puff_t <- ifelse(is.na(iM0$TTLtime[1]), puff_default, iM0$TTLtime[1])
#
#   key <- paste(
#     typ,
#     paste0("subc=", iM0$assigned_subclass[1]),
#     paste0("nSweep=", n_sweeps),
#     paste0("nPuff=", length(chunk_info$puff_pos)),
#     paste0("sweeps=", trim1(iM0$wrapped_sweeps[1])),
#     paste0("puffTime=", puff_t),
#     paste0("bath=", safe_chr(iM0$Bath_Drug[1])),
#     paste0("puffdrug=", safe_chr(iM0$Puff_drug[1])),
#     paste0("exp=", safe_chr(iM0$Exp[1])),
#     sep = "|"
#   )
#
#   gsub("[^A-Za-z0-9|=_\\-]+", "_", key)
# }
#
# # Build aligned raw trace + segment_bounds
# # - Uses first `tpRMP_puff` sweep as the "zero reference" for concatenation
# # - Additionally shifts ALL time by puffTime so puff onset lands near 0 for the puff sweep
# build_raw_trace_rmp <- function(sweep_objs, dur, chunk_names, puffTime) {
#   # sweep_objs: named list; each element is list(df=..., seg=...)
#   #names(sweep_objs) <- chunk_names
#   #names(dur) <- chunk_names
#
#   sweep_dfs <- lapply(sweep_objs, `[[`, "df")
#
#   cum_shift <- c(0, cumsum(as.numeric(dur)))[seq_along(sweep_dfs)]
#
#   first_puff_idx <- which(chunk_names == "tpRMP_puff")[1]
#   if (is.na(first_puff_idx)) first_puff_idx <- 1L
#   zero_offset <- cum_shift[first_puff_idx]
#
#   segment_bounds <- data.frame(
#     sweep = seq_along(sweep_dfs),
#     chunk = chunk_names,
#     start_time = NA_real_,
#     end_time = NA_real_,
#     puff_onset_time = NA_real_,
#     stringsAsFactors = FALSE
#   )
#
#   out_list <- vector("list", length(sweep_dfs))
#
#   for (i in seq_along(sweep_dfs)) {
#     d <- sweep_dfs[[i]]
#
#     # global concat in seconds
#     t_pre <- d$epoch_t + cum_shift[i] - zero_offset
#     t <- t_pre - puffTime
#
#     out_list[[i]] <- data.frame(
#       time = t,
#       tpRMP = d$tpRMP,
#       pulse_state = d$pulse_state %||% NA_character_,
#       chunk = chunk_names[i],
#       sweep = i
#     )
#
#     segment_bounds$start_time[i] <- min(t, na.rm = TRUE)
#     segment_bounds$end_time[i]   <- max(t, na.rm = TRUE)
#
#     if (chunk_names[i] == "tpRMP_puff" && nrow(d) > 0) {
#       segment_bounds$puff_onset_time[i] <- 0
#     }
#   }
#
#   raw <- do.call(rbind, out_list)
#   list(raw_trace = raw, segment_bounds = segment_bounds)
# }
#
# analyze_pulse_off_rmp <- function(tst, sr = NULL) {
#   # tst is your rmp_puff_output table with columns time, tpRMP, pulse_state
#   stopifnot(all(c("time","tpRMP","pulse_state") %in% names(tst)))
#
#   d_off <- tst[tst$pulse_state == "pulse_off" & is.finite(tst$tpRMP), , drop = FALSE]
#   if (nrow(d_off) == 0) {
#     return(data.frame(
#       n = 0, mean_mV = NA_real_, sd_mV = NA_real_, median_mV = NA_real_,
#       p05 = NA_real_, p95 = NA_real_
#     ))
#   }
#
#   data.frame(
#     n = nrow(d_off),
#     mean_mV = mean(d_off$tpRMP),
#     sd_mV = stats::sd(d_off$tpRMP),
#     median_mV = stats::median(d_off$tpRMP),
#     p05 = stats::quantile(d_off$tpRMP, 0.05, names = FALSE),
#     p95 = stats::quantile(d_off$tpRMP, 0.95, names = FALSE)
#   )
# }
#
# analyze_pulse_off_sections <- function(tst,
#                                        time_col = "time",
#                                        y_col = "tpRMP",
#                                        state_col = "pulse_state",
#                                        on_label = "pulse_on",
#                                        off_label = "pulse_off",
#                                        max_gap_s = NULL,
#                                        sampling_rate = NULL,
#                                        min_duration_s = 0) {
#
#   stopifnot(is.data.frame(tst))
#   stopifnot(all(c(time_col, y_col, state_col) %in% names(tst)))
#
#   d <- tst
#   d <- d[is.finite(d[[time_col]]) & is.finite(d[[y_col]]), , drop = FALSE]
#   d <- d[order(d[[time_col]]), , drop = FALSE]
#
#   d_off <- d[d[[state_col]] == off_label, , drop = FALSE]
#   if (nrow(d_off) == 0) {
#     return(list(
#       sections = data.frame(),
#       overall = data.frame(
#         n = 0, mean_mV = NA_real_, sd_mV = NA_real_, median_mV = NA_real_,
#         p05 = NA_real_, p95 = NA_real_
#       )
#     ))
#   }
#
#   # Determine a default max gap tolerance if not provided
#   # (This helps prevent tiny missing samples from splitting a section.)
#   if (is.null(max_gap_s)) {
#     if (!is.null(sampling_rate) && is.numeric(sampling_rate) && sampling_rate > 0) {
#       max_gap_s <- 1.5 / sampling_rate   # tolerate ~1 missed sample
#     } else {
#       # fall back to an empirical step estimate
#       dt_med <- stats::median(diff(d_off[[time_col]]), na.rm = TRUE)
#       max_gap_s <- 1.5 * dt_med
#     }
#   }
#
#   # Identify contiguous "off" sections by time gaps
#   dt <- c(Inf, diff(d_off[[time_col]]))
#   new_section <- dt > max_gap_s
#   section_id <- cumsum(new_section)
#
#   d_off$off_section_id <- section_id
#
#   # Summarize each off section
#   sections <- d_off |>
#     dplyr::group_by(off_section_id) |>
#     dplyr::summarise(
#       t_start = min(.data[[time_col]], na.rm = TRUE),
#       t_end   = max(.data[[time_col]], na.rm = TRUE),
#       duration_s = t_end - t_start,
#       n = dplyr::n(),
#       mean_mV = mean(.data[[y_col]], na.rm = TRUE),
#       sd_mV = stats::sd(.data[[y_col]], na.rm = TRUE),
#       median_mV = stats::median(.data[[y_col]], na.rm = TRUE),
#       p05 = stats::quantile(.data[[y_col]], 0.05, names = FALSE, na.rm = TRUE),
#       p95 = stats::quantile(.data[[y_col]], 0.95, names = FALSE, na.rm = TRUE),
#       .groups = "drop"
#     ) |>
#     dplyr::filter(duration_s >= min_duration_s) |>
#     dplyr::mutate(off_section_id = as.integer(off_section_id)) |>
#     dplyr::arrange(t_start)
#
#   # Overall pooled (optional but often useful)
#   overall <- data.frame(
#     n = nrow(d_off),
#     mean_mV = mean(d_off[[y_col]], na.rm = TRUE),
#     sd_mV = stats::sd(d_off[[y_col]], na.rm = TRUE),
#     median_mV = stats::median(d_off[[y_col]], na.rm = TRUE),
#     p05 = stats::quantile(d_off[[y_col]], 0.05, names = FALSE, na.rm = TRUE),
#     p95 = stats::quantile(d_off[[y_col]], 0.95, names = FALSE, na.rm = TRUE)
#   )
#
#   list(sections = sections, overall = overall)
# }
