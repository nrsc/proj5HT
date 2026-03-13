#' Plot concatenated ASP protocols with a trace + heat strip
#'
#' Concatenates all protocol entries in `asp$by_protocol` into one continuous
#' **global time** axis separated by a fixed gap (`gap_s`). Produces a two-panel
#' plot: (1) the response trace (points colored by `instRate`) and (2) a one-row
#' heat strip of `instRate` binned by `heat_dt`, with green TTL markers placed
#' at each protocol’s TTL onset (converted into global time).
#'
#' Protocols can be ordered either by numeric sweep indices parsed from
#' `pro_row$wrapped_sweeps` or by `protocol_key`. Gaps between protocols may be
#' indicated by dotted vertical lines at each protocol start.
#'
#' @param asp A nested list returned by `build_asp()`. Must contain `asp$by_protocol`
#'   (named list) and `asp$md$cell_name`. Each `asp$by_protocol[[key]]` should include
#'   `spike_puff_output` with columns `time` and `instRate`, and optionally
#'   `percent_change`, `blMean`, `segment_bounds`, and `pro_row` metadata used for
#'   ordering/labels (e.g. `wrapped_sweeps`, `stimulus_description`, `Exp`,
#'   `assigned_subclass`).
#' @param gap_s Numeric. Gap duration (seconds) inserted between consecutive protocol
#'   segments in global time.
#' @param order_by Character. Ordering rule for protocols:
#'   \describe{
#'     \item{`"wrapped_sweeps"`}{Order by numeric sweep ids parsed from `pro_row$wrapped_sweeps`
#'       (min sweep, then max sweep, then type/key).}
#'     \item{`"protocol_key"`}{Order by protocol type then key.}
#'   }
#' @param bl_window Numeric length-2 vector. Baseline window in seconds (e.g. `c(-5,0)`),
#'   used when `percent_change` is missing and/or `blMean` is invalid.
#' @param y_mode Character. Y-axis transform for the trace panel:
#'   \describe{
#'     \item{`"percent"`}{plot `% of baseline` (`percent_change`), with reference line at 100}
#'     \item{`"log10_percent"`}{plot `log10(percent_change)` (drops non-positive values), reference at `log10(100)`}
#'     \item{`"signed_log10_delta"`}{plot `sign(%-100) * log10(|%-100| + 1)`, reference at 0}
#'     \item{`"log2_fc"`}{plot `log2(pmax(instRate, min_rate_hz) / baseline)`, reference at 0}
#'   }
#' @param y_arrow_pad Optional numeric. Vertical padding added above each protocol’s
#'   segment maximum (in `y_plot` units) to place TTL markers. If `NULL`, a default
#'   is chosen based on `y_mode`.
#' @param y_lim Numeric length-2 vector. Visible y-range for the trace panel via
#'   `coord_cartesian(ylim = y_lim)`.
#' @param alpha_trace Numeric in [0,1]. Alpha for trace points in the main panel.
#' @param linewidth Numeric. Line width for the faint black backbone trace.
#' @param show_gap_lines Logical. If `TRUE`, draw dotted vertical lines at each
#'   protocol segment start (`df_starts$x0`) in both panels.
#' @param gap_line_alpha Numeric in [0,1]. Alpha for the dotted gap/protocol-start lines.
#' @param min_rate_hz Numeric. Minimum firing rate (Hz) used as a floor when
#'   `y_mode = "log2_fc"` to prevent `-Inf` and preserve the trace shape through
#'   silent periods.
#' @param heat_dt Numeric. Bin width (seconds) for the heat strip. InstRate is
#'   aggregated within bins and filled to create a continuous strip.
#' @param heat_method Character. Aggregation for instRate within each bin: `"mean"` or
#'   `"median"`.
#' @param save_png Logical. If `TRUE`, save the combined plot to disk via `ggsave()`.
#' @param rookery_path Character. Base directory used when constructing a default
#'   `png_path` (typically a per-cell folder under `"rookery"`).
#' @param png_path Character. Output file path for the PNG. If `NULL` and
#'   `save_png = TRUE`, defaults to
#'   `file.path(rookery_path, asp$md$cell_name, paste0(asp$md$cell_name, "-asp-concat-heat.png"))`.
#'   If a directory path is provided (no `.png` suffix), the file name is appended.
#' @param png_w Numeric. Width (in inches) passed to `ggsave()`.
#' @param png_h Numeric. Height (in inches) passed to `ggsave()`.
#' @param dpi Integer. DPI passed to `ggsave()`.
#' @param verbose Logical. If `TRUE`, emit per-protocol TTL diagnostics and a TTL
#'   sanity check summary.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{`df_concat`}{Concatenated trace data for all plotted protocols, including
#'     `time_global`, `protocol_key`, `protocol_index`, and `y_plot`.}
#'   \item{`df_arrows`}{TTL marker table in global time with per-protocol y placement.}
#'   \item{`df_starts`}{Protocol segment start/end table (`x0`, `x1`) in global time.}
#'   \item{`df_heat`}{Heat-strip binned data (global time bins and instRate). Bins
#'     outside any protocol segment are set to `NA`.}
#'   \item{`ttl_check`}{Per-protocol TTL sanity check summary (source, counts, range check).}
#'   \item{`plot`}{The combined patchwork plot object.}
#' }
#'
#' @details
#' **Global time construction:** each protocol’s local time is shifted by `seg_start`
#' such that the segment begins at the current `offset`. The offset is then advanced
#' to the end of the segment plus `gap_s`.
#'
#' **Percent change:** if `percent_change` is absent, it is computed as
#' `(instRate / baseline) * 100` using `obj$blMean` when valid, otherwise from
#' `bl_window`.
#'
#' **TTL placement:** TTL times are sourced in this priority order:
#' \enumerate{
#'   \item `segment_bounds$ttl_onset_time` for TTL-related chunks when present and finite
#'   \item `segment_bounds$start_time + obj$ttlTime` as a fallback when onset times are absent
#'   \item local `t = 0` when no usable segment bounds are available
#' }
#' TTL times are then snapped to the nearest available sample time in the local trace
#' and shifted into global time for plotting.
#'
#' **Heat strip:** bins are computed across the full global span, filled in both
#' directions, and then blanked (`NA`) outside any protocol segment.
#'
#' @examples
#' \dontrun{
#' # Basic concatenated view ordered by sweeps
#' plot_asp_concatenated(asp, order_by = "wrapped_sweeps", y_mode = "log2_fc")
#'
#' # Tight y-range, larger gaps, save to rookery
#' plot_asp_concatenated(
#'   asp,
#'   gap_s = 3,
#'   y_mode = "signed_log10_delta",
#'   y_lim = c(-1.5, 1.5),
#'   save_png = TRUE
#' )
#'
#' # Use median heat bins and hide gap lines
#' plot_asp_concatenated(
#'   asp,
#'   heat_method = "median",
#'   show_gap_lines = FALSE
#' )
#'
#' order_by = "wrapped_sweeps"
#' y_mode = "log2_fc"
#' bl_window = c(-5, 0)
#' y_arrow_pad = NULL
#' y_lim = c(-1,1)
#' alpha_trace = 0.9
#' linewidth = 0.35
#' show_gap_lines = TRUE
#' min_rate_hz = 0.05
#' gap_line_alpha = 0.35
#' heat_dt = 0.10
#' heat_method = c("mean", "median")
#' verbose = TRUE
#'
#' }
#' @export
plot_asp_concatenated <- function(asp,
                                  gap_s = 2,
                                  order_by = c("wrapped_sweeps", "protocol_key"),
                                  bl_window = c(-5, 0),
                                  y_mode = c("percent", "log10_percent", "signed_log10_delta", "log2_fc"),
                                  y_arrow_pad = NULL,
                                  y_lim = NULL,
                                  x_lim = NULL,
                                  x_lim_by_protocol = NULL,
                                  alpha_trace = 0.9,
                                  linewidth = 0.35,
                                  show_gap_lines = TRUE,
                                  gap_line_alpha = 0.35,
                                  min_rate_hz = 0.05,
                                  heat_dt = 0.10,
                                  heat_method = c("mean", "median"),
                                  save_png = FALSE,
                                  rookery_path = "rookery",
                                  png_path = NULL,
                                  png_w = 14,
                                  png_h = 5,
                                  dpi = 300,
                                  verbose = TRUE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
  })

  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  resolve_xlim_local <- function(key, x_lim_by_protocol) {
    if (is.null(x_lim_by_protocol)) return(NULL)

    # Case A: a single numeric length-2 vector applies to ALL protocols
    if (is.numeric(x_lim_by_protocol) && length(x_lim_by_protocol) == 2) {
      if (all(is.finite(x_lim_by_protocol))) return(as.numeric(x_lim_by_protocol))
      return(NULL)
    }

    # Case B: a named list, e.g. list("keyA"=c(-5,30), "keyB"=c(0,60))
    if (is.list(x_lim_by_protocol) && !is.data.frame(x_lim_by_protocol)) {
      if (!is.null(names(x_lim_by_protocol)) && key %in% names(x_lim_by_protocol)) {
        xl <- x_lim_by_protocol[[key]]
        if (is.numeric(xl) && length(xl) == 2 && all(is.finite(xl))) return(as.numeric(xl))
      }
      return(NULL)
    }

    # Case C: a data.frame with columns protocol_key, xmin, xmax
    if (is.data.frame(x_lim_by_protocol)) {
      if (all(c("protocol_key", "xmin", "xmax") %in% names(x_lim_by_protocol))) {
        row <- x_lim_by_protocol[x_lim_by_protocol$protocol_key == key, , drop = FALSE]
        if (nrow(row) >= 1) {
          xl <- c(row$xmin[1], row$xmax[1])
          if (is.numeric(xl) && length(xl) == 2 && all(is.finite(xl))) return(as.numeric(xl))
        }
      }
      return(NULL)
    }

    NULL
  }


  parse_exp_from_key <- function(key) {
    key <- as.character(key)[1]  # force scalar
    m <- sub(".*\\|exp=([^|]*)(\\||$).*", "\\1", key)
    m <- trimws(m)

    if (is.na(m) || m == "" || identical(m, key)) {
      "standard puff"
    } else if (grepl("wash-in", m, ignore.case = TRUE)) {
      "wash-in"
    } else if (grepl("wash-out|wout", m, ignore.case = TRUE)) {
      "wash-out"
    } else {
      m
    }
  }


  rank_label <- function(obj) {
    ex <- obj$pro_row$Rank

    if (length(ex) == 0) return("")

    ex <- trimws(as.character(ex)[1])

    if (is.na(ex) || ex == "") {
      return("")
    }

    if (ex == "1" || ex == "blSpike") {
      return("Representative Baseline")
    }

    if (ex == "testSpike") {
      return("Test Trace")
    }

    ex
  }

  exp_label <- function(obj) {
    ex <- obj$pro_row$Exp
    ex <- if (length(ex) == 0) NA_character_ else as.character(ex)[1]
    ex <- trimws(ex)
    if (is.na(ex) || ex == "") "none" else ex
  }

  bath_label <- function(obj) {
    ex <- obj$pro_row$Bath_Drug
    ex <- if (length(ex) == 0) NA_character_ else as.character(ex)[1]
    ex <- trimws(ex)
    if (is.na(ex) || ex == "") "none" else ex
  }

  order_by <- match.arg(order_by)
  y_mode   <- match.arg(y_mode)
  heat_method <- match.arg(heat_method)

  if (is.null(asp) || is.null(asp$by_protocol) || length(asp$by_protocol) == 0) {
    stop("asp$by_protocol is empty.")
  }

  # -----------------------------
  # meta table for ordering (by wrapped_sweeps numeric values)
  # -----------------------------
  parse_sweeps_vec <- function(w) {
    w <- as.character(w)
    # extract all integers appearing in the string
    nums <- suppressWarnings(as.numeric(unlist(regmatches(w, gregexpr("[0-9]+", w)))))
    nums <- nums[is.finite(nums)]
    nums
  }

  meta <- dplyr::bind_rows(lapply(names(asp$by_protocol), function(key) {
    obj <- asp$by_protocol[[key]]
    pr  <- obj$pro_row

    stim <- as.character(pr$stimulus_description)
    proto_type <- dplyr::case_when(
      grepl("TTLstack", stim) ~ "TTLstack",
      grepl("spikeTTL", stim) ~ "spikeTTL",
      TRUE ~ "other"
    )

    sweeps_vec <- parse_sweeps_vec(pr$wrapped_sweeps)
    min_sweep <- if (length(sweeps_vec)) min(sweeps_vec) else Inf
    max_sweep <- if (length(sweeps_vec)) max(sweeps_vec) else Inf

    data.frame(
      protocol_key = key,
      protocol_type = proto_type,
      stimulus_description = stim,
      wrapped_sweeps = as.character(pr$wrapped_sweeps),
      min_sweep = min_sweep,
      max_sweep = max_sweep,
      stringsAsFactors = FALSE
    )
  }))

  if (order_by == "wrapped_sweeps") {
    meta <- meta %>%
      dplyr::arrange(.data$min_sweep, .data$max_sweep, .data$protocol_type, .data$protocol_key)
  } else {
    meta <- meta %>%
      dplyr::arrange(.data$protocol_type, .data$protocol_key)
  }

  # -----------------------------
  # title subtitle (pred + assigned)
  # -----------------------------
  assigned_subs <- unique(na.omit(vapply(
    asp$by_protocol,
    function(obj) as.character(obj$pro_row$assigned_subclass),
    character(1)
  )))

  assigned_sub_txt <- if (length(assigned_subs) == 0) {
    "assigned=NA"
  } else if (length(assigned_subs) == 1) {
    paste0("assigned=", assigned_subs)
  } else {
    paste0("assigned={", paste(assigned_subs, collapse = ","), "}")
  }

  pred_txt <- if (!is.null(asp$md$predicted_subclass) && !is.na(asp$md$predicted_subclass)) {
    paste0("predicted=", asp$md$predicted_subclass)
  } else {
    "predicted=NA"
  }

  blockers_txt <- if (!is.null(asp$md$Fast_blockers) && !is.na(asp$md$Fast_blockers)) {
    paste0("blockers=", asp$md$Fast_blockers)
  } else {
    "blockers=NA"
  }

  subtitle_txt <- paste(assigned_sub_txt, blockers_txt, sep = " | ")

  # -----------------------------
  # y arrow pad default per mode
  # -----------------------------
  if (is.null(y_arrow_pad)) {
    y_arrow_pad <- switch(
      y_mode,
      percent = 5,
      log10_percent = 0.06,          # pad expressed in raw % units (we place dots in transformed y later)
      signed_log10_delta = 0.20,
      log2_fc = 0.25
    )
  }

  # -----------------------------
  # concatenate traces + TTL markers + starts
  # -----------------------------
  offset <- 0
  out_list <- list()
  arrow_list <- list()
  starts_list <- list()
  ttl_check <- list()
  n_fallback_t0 <- 0

  for (ii in seq_len(nrow(meta))) {

    key <- meta$protocol_key[ii]
    obj <- asp$by_protocol[[key]]
    pr  <- obj$pro_row

    d <- obj$spike_puff_output
    if (is.null(d) || !"time" %in% names(d)) next

    # require time + instRate at minimum
    if (!"instRate" %in% names(d)) {
      if (verbose) message("Skipping ", key, ": missing instRate in spike_puff_output.")
      next
    }

    # compute percent_change if missing
    if (!"percent_change" %in% names(d)) {
      if (verbose) message("Protocol ", key, " missing percent_change; computing from instRate.")

      blMean <- obj$blMean
      if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
        bl_keep <- d$time > bl_window[1] & d$time <= bl_window[2]
        blMean <- mean(d$instRate[bl_keep], na.rm = TRUE)
      }

      if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
        if (verbose) message("Skipping ", key, ": cannot compute percent_change (bad blMean).")
        next
      }

      d$percent_change <- (d$instRate / blMean) * 100
    }

    # Now safe to filter
    d <- d %>%
      dplyr::filter(is.finite(time), is.finite(instRate), is.finite(percent_change)) %>%
      dplyr::arrange(time)

    # -----------------------------
    # OPTIONAL: per-protocol local x-limits (crop BEFORE concatenation)
    # -----------------------------
    xlim_local <- resolve_xlim_local(key, x_lim_by_protocol)
    if (!is.null(xlim_local)) {
      d <- dplyr::filter(d, .data$time >= xlim_local[1], .data$time <= xlim_local[2])
      d <- dplyr::arrange(d, .data$time)
    }
    if (nrow(d) == 0) next

    # global time placement
    tmin <- min(d$time, na.rm = TRUE)
    tmax <- max(d$time, na.rm = TRUE)

    seg_start <- offset - tmin
    d$time_global <- d$time + seg_start
    d$protocol_key <- key
    d$protocol_index <- ii
    d$protocol_type <- meta$protocol_type[ii]
    d$stimulus_description <- meta$stimulus_description[ii]



    # compute y_plot
    if (y_mode == "percent") {
      d$y_plot <- d$percent_change
      y_lab <- "% of baseline instRate"
      y_hline <- 100
    } else if (y_mode == "log10_percent") {
      d <- d %>% dplyr::filter(percent_change > 0)
      d$y_plot <- log10(d$percent_change)
      y_lab <- "log10(% of baseline)"
      y_hline <- log10(100)
    } else if (y_mode == "signed_log10_delta") {
      delta <- d$percent_change - 100
      d$y_plot <- sign(delta) * log10(abs(delta) + 1)
      y_lab <- "sign(%-100) * log10(|%-100|+1)"
      y_hline <- 0
    } else if (y_mode == "log2_fc") {

      # Use a pseudocount for instRate so log2() never sees 0.
      # This preserves the "line/shape" through silent periods without shooting to -Inf.
      blMean <- obj$blMean
      if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
        bl_keep <- d$time > bl_window[1] & d$time <= bl_window[2]
        blMean <- mean(d$instRate[bl_keep], na.rm = TRUE)
      }
      if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
        if (verbose) message("Skipping ", key, ": cannot compute log2_fc (bad blMean).")
        next
      }

      d$y_plot <- log2(pmax(d$instRate, min_rate_hz) / blMean)
      y_lab <- "log2(instRate / baseline)"
      y_hline <- 0
    }


    if (nrow(d) == 0) next
    out_list[[key]] <- d

    # protocol start (dotted line)
    starts_list[[key]] <- data.frame(
      protocol_key = key,
      protocol_index = ii,
      x0 = seg_start + tmin,  # start of protocol in global time
      x1 = seg_start + tmax,  # end of protocol in global time
      stringsAsFactors = FALSE
    )

    # TTL times from segment_bounds if available, else fallback to local t=0
    ttl_times_local <- NULL
    ttl_source <- "fallback_t0"

    if (!is.null(obj$segment_bounds) &&
        all(c("chunk", "start_time") %in% names(obj$segment_bounds))) {

      sb <- obj$segment_bounds
      ttl_rows <- which(sb$chunk == "spikeTTL" | grepl("^TTLstack", sb$chunk))

      if (length(ttl_rows) > 0 &&
          "ttl_onset_time" %in% names(sb) &&
          all(is.finite(sb$ttl_onset_time[ttl_rows]))) {

        ttl_times_local <- sb$ttl_onset_time[ttl_rows]
        ttl_source <- "segment_bounds_ttl_onset"

      } else if (length(ttl_rows) > 0 &&
                 all(is.finite(sb$start_time[ttl_rows]))) {

        # fallback: infer TTL onset from sweep start + ttlTime
        ttl_times_local <- sb$start_time[ttl_rows] + obj$ttlTime
        ttl_source <- "segment_bounds_start_plus_ttlTime"
      }

    }

    if (is.null(ttl_times_local) || length(ttl_times_local) == 0) {
      ttl_times_local <- 0
      n_fallback_t0 <- n_fallback_t0 + 1
    }

    # Snap TTL times to nearest available trace timepoint
    ttl_times_local <- vapply(ttl_times_local, function(tt) {
      d$time[which.min(abs(d$time - tt))]
    }, numeric(1))


    if (verbose) {
      message("TTL source=", ttl_source, " | ttl_times_local=",
              paste(round(ttl_times_local, 3), collapse = ","))
    }


    ttl_times_global <- ttl_times_local + seg_start
    if (verbose) {
      message(
        "TTL debug | ", key,
        " | seg_start=", round(seg_start, 3),
        " | tmin=", round(tmin, 3),
        " | tmax=", round(tmax, 3),
        " | ttl_local=", paste(round(ttl_times_local, 3), collapse = ","),
        " | ttl_global=", paste(round(ttl_times_global, 3), collapse = ",")
      )
    }

    # TTL sanity check: ensure within local trace time span
    ttl_ok <- ttl_times_local >= min(d$time, na.rm = TRUE) - 1e-6 &
      ttl_times_local <= max(d$time, na.rm = TRUE) + 1e-6

    ttl_check[[key]] <- data.frame(
      protocol_key = key,
      protocol_index = ii,
      ttl_source = ttl_source,
      n_ttl = length(ttl_times_local),
      ttl_ok_all = all(ttl_ok),
      stringsAsFactors = FALSE
    )

    # place TTL dots in y_plot space
    y_seg_max <- max(d$y_plot, na.rm = TRUE)
    y_dot <- y_seg_max + y_arrow_pad

    arrow_list[[key]] <- data.frame(
      protocol_key   = key,
      protocol_index = ii,
      time_global    = ttl_times_global,   # <-- correct
      y_dot          = y_dot,
      stringsAsFactors = FALSE
    )

    # advance offset (end of segment + gap)
    if (ii < nrow(meta)) {
      offset <- seg_start + tmax + gap_s
    } else {
      offset <- seg_start + tmax
    }

  }

  df_starts <- dplyr::bind_rows(starts_list) %>%
    dplyr::arrange(.data$x0)



  # -----------------------------
  # Protocol EXP labels (AFTER df_starts exists)
  # -----------------------------
  get_exp_label <- function(obj) {
    ex <- obj$pro_row$Exp
    ex <- if (length(ex) == 0) NA_character_ else as.character(ex)[1]
    ex <- trimws(ex)
    if (is.na(ex) || ex == "") "baseline puff" else ex
  }

  df_labels <- dplyr::bind_rows(lapply(seq_len(nrow(meta)), function(ii) {
    key <- meta$protocol_key[ii]
    obj <- asp$by_protocol[[key]]
    if (is.null(obj)) return(NULL)

    st <- df_starts %>% dplyr::filter(.data$protocol_key == key)
    if (nrow(st) == 0) return(NULL)

    data.frame(
      protocol_key = key,
      protocol_index = ii,
      x = (st$x0[1] + st$x1[1]) / 2,
      exp_label = paste(get_exp_label(obj), rank = rank_label(obj), sep = " | "),
      stringsAsFactors = FALSE
    )
  }))


  # df_gaps <- df_starts %>%
  #   dplyr::mutate(
  #     gap_start = .data$x1,
  #     gap_end   = dplyr::lead(.data$x0)
  #   ) %>%
  #   dplyr::filter(is.finite(.data$gap_start), is.finite(.data$gap_end), .data$gap_end > .data$gap_start)


  df_concat <- bind_rows(out_list)

  #range(df_concat$instRate)

  df_arrows <- bind_rows(arrow_list)
  df_ttlchk <- bind_rows(ttl_check)

  #df_arrows

  if (nrow(df_concat) == 0) stop("No usable traces found after processing protocols.")

  # shared x limits + instRate colour limits (define BEFORE plotting)
  xlim_global <- range(df_concat$time_global, na.rm = TRUE)
  xlim_use <- if (!is.null(x_lim)) x_lim else xlim_global
  rate_lim <- range(df_concat$instRate, na.rm = TRUE)
  max_rate <- max(df_concat$instRate, na.rm = TRUE)

  if (!is.null(y_lim)) {
    ylim_use <- y_lim
  } else {
    yr <- range(d$y_plot, na.rm = TRUE)
    ylim_use <- c(max(yr[1], -1.5), yr[2])
  }

  if (isTRUE(verbose)) {
    bad <- df_ttlchk %>% dplyr::filter(!ttl_ok_all)
    if (nrow(bad) > 0) {
      message("TTL sanity check: some protocols have TTL times outside trace range:")
      print(bad)
    } else {
      message("TTL sanity check: OK (TTL times within trace windows).")
    }
    if (n_fallback_t0 > 0) {
      message("TTL timing fallback used for ", n_fallback_t0,
              " protocol(s) (no segment_bounds; placing TTL at local t=0).")
    }
  }

  p_trace_title = paste(asp$md$cell_name, "concatenated spike protocols",  sep = " | ")
  # -----------------------------
  # Main trace plot
  # -----------------------------
  p_trace <- ggplot(df_concat, aes(x = time_global, y = y_plot)) +
    geom_hline(yintercept = y_hline, linewidth = 0.3) +

    # optional: keep a faint black backbone so the shape is readable
    geom_line(colour = "black", alpha = 0.25, linewidth = linewidth) +

    # main: points colored by instRate
    geom_point(aes(colour = instRate),
               alpha = alpha_trace,
               size = 0.30) +

    scale_colour_gradient(
      low = "navy",
      high = "orange",
      limits = c(0, max_rate),
      oob = scales::squish,
      na.value = "navy") +


    guides(colour = "none") +

    labs(
      title = p_trace_title,
      subtitle = subtitle_txt,
      x = NULL,
      y = y_lab
    ) +
    theme_classic() +
    theme(plot.title.position = "plot")


  p_trace <- p_trace +
    scale_x_continuous(limits = xlim_global, expand = c(0, 0)) +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    theme(plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

  # Put labels near the top of the trace panel
    if (!is.null(y_lim) && length(y_lim) == 2 && all(is.finite(y_lim))) {
      y_lab_pos <- y_lim[2] - 0.05 * diff(y_lim)
    } else {
      yr <- range(df_concat$y_plot, na.rm = TRUE)
      y_lab_pos <- yr[2] - 0.05 * diff(yr)
    }

  if (nrow(df_labels) > 0) {
    p_trace <- p_trace +
      geom_label(
        data = df_labels,
        aes(x = x, y = y_lab_pos, label = exp_label),
        inherit.aes = FALSE,
        size = 2.6,
        label.size = 0.2,
        label.padding = grid::unit(0.10, "lines"),
        label.r = grid::unit(0.12, "lines"),
        alpha = 0.85
      )
  }





  if (isTRUE(show_gap_lines) && nrow(df_starts) > 0) {
    p_trace <- p_trace +
      geom_vline(
        data = df_starts,
        aes(xintercept = x0),
        linetype = "dotted",
        linewidth = 0.4,
        alpha = gap_line_alpha#,
        #inherit.aes = FALSE
      )
  }

  # -----------------------------
  # Heatmap strip: frequency (instRate) over concatenated time
  # -----------------------------
  # Explicitly ensure instRate exists here (it should, from df_concat)
  if (!"instRate" %in% names(df_concat)) {
    stop("Internal error: df_concat lacks instRate, cannot build heatmap.")
  }

  xlim_global <- range(df_concat$time_global, na.rm = TRUE)

  # 1) compute bin summaries from the (possibly sparse) data
  df_heat_raw <- df_concat %>%
    dplyr::mutate(time_bin_start = floor(.data$time_global / heat_dt) * heat_dt) %>%
    dplyr::group_by(.data$time_bin_start) %>%
    dplyr::summarise(
      instRate = if (heat_method == "mean") mean(.data$instRate, na.rm = TRUE) else median(.data$instRate, na.rm = TRUE),
      .groups = "drop"
    )

  # # 2) make a full bin grid across the entire global range
  # df_heat_raw <- df_concat %>%
  #   dplyr::mutate(time_bin_start = floor(.data$time_global / heat_dt) * heat_dt) %>%
  #   dplyr::group_by(.data$time_bin_start) %>%
  #   dplyr::summarise(
  #     instRate = if (heat_method == "mean") mean(.data$instRate, na.rm = TRUE) else median(.data$instRate, na.rm = TRUE),
  #     .groups = "drop"
  #   )

  bin_seq <- seq(
    from = floor(xlim_global[1] / heat_dt) * heat_dt,
    to   = floor(xlim_global[2] / heat_dt) * heat_dt,
    by   = heat_dt
  )

  df_heat <- tibble::tibble(time_bin_start = bin_seq) %>%
    dplyr::left_join(df_heat_raw, by = "time_bin_start") %>%
    dplyr::mutate(
      time_bin = .data$time_bin_start + heat_dt/2,
      y = 1
    )

  df_heat$instRate[is.na(df_heat$instRate)] = 0

  # mark bins inside any protocol segment (no eps trimming)
  in_any_segment <- rep(FALSE, nrow(df_heat))
  for (j in seq_len(nrow(df_starts))) {
    in_any_segment <- in_any_segment |
      (df_heat$time_bin >= df_starts$x0[j] & df_heat$time_bin <= df_starts$x1[j])
  }
  df_heat$instRate[!in_any_segment] <- NA_real_


  p_heat <- ggplot(df_heat, aes(x = time_bin, y = y, fill = instRate)) +
    geom_tile(height = 1, width = heat_dt) +
    scale_x_continuous(limits = xlim_global, expand = c(0, 0)) +
    coord_cartesian(xlim = xlim_use, ylim = c(0.5,1.8), expand = FALSE) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      plot.margin = margin(t = 0, r = 5, b = 5, l = 5),
      legend.position = "right"
    ) +
    labs(
      x = "Global time (s) (protocols concatenated with gaps)",
      y = NULL,
      fill = "instRate (Hz)"
    )


  if (isTRUE(show_gap_lines) && nrow(df_starts) > 0) {
    p_heat <- p_heat +
      geom_vline(
        data = df_starts,
        aes(xintercept = x0),
        linetype = "dotted",
        linewidth = 0.4,
        alpha = gap_line_alpha#,
        #inherit.aes = FALSE
      )
  }

  max_rate <- max(df_concat$instRate, na.rm = TRUE)

  p_heat <- p_heat +
    scale_fill_gradient(
      low = "navy",          # darker than "blue"
      high = "orange",
      limits = c(0, max_rate),
      oob = scales::squish,
      na.value = "white"    # for NA (gaps outside protocol segments)
    )

  if (nrow(df_arrows) > 0) {
    p_heat <- p_heat +
      geom_point(
        data = df_arrows,
        aes(x = time_global, y = 1.62),
        inherit.aes = FALSE,
        shape = 6,
        size = 2,
        colour = "green4"
      )
  }



  # -----------------------------
  # Combine
  # -----------------------------
  p <- p_trace / p_heat + patchwork::plot_layout(heights = c(3, 1))

  print(p)

  if (isTRUE(save_png)) {
    if (is.null(png_path)) {
      png_path <- file.path(rookery_path, asp$md$cell_name, paste0(asp$md$cell_name, "-asp-concat-heat.png"))
    }
    if(!grepl(".png", png_path)){
      png_path <- file.path(png_path, paste0(asp$md$cell_name, "-asp-concat-heat.png"))
    }

    ggsave(filename = png_path, plot = p, width = png_w, height = png_h, dpi = dpi)
    if (isTRUE(verbose)) message("Saved: ", png_path)
  }

  invisible(list(
    df_concat = df_concat,
    df_arrows = df_arrows,
    df_starts = df_starts,
    df_heat = df_heat,
    ttl_check = df_ttlchk,
    plot = p
  ))
}
