#' Plot a single protocol section from an ASP object
#'
#' Builds a two-panel plot for **one** protocol entry inside `asp$by_protocol`:
#' (1) a trace of the response over time (with points colored by instantaneous firing rate),
#' and (2) a one-row heat strip (binned by `heat_dt`) showing instRate over time with a
#' green down-arrow at the TTL onset.
#'
#' This function assumes the selected protocol's `spike_puff_output$time` is already
#' in **local time** with TTL onset near **0 s** (i.e., do *not* use global/concatenated
#' segment start times to place the TTL marker).
#'
#' @param asp A nested list returned by `build_asp()`. Must contain `asp$by_protocol`
#'   (named list) and `asp$md$cell_name`. Each entry in `asp$by_protocol[[key]]` must
#'   include `spike_puff_output` with at least `time` and `instRate`, and optionally
#'   `percent_change`, `blMean`, and `segment_bounds`.
#' @param protocol_key Character. Exact name of the protocol entry to plot (must match
#'   one of `names(asp$by_protocol)`). If provided, takes precedence over `key_match`.
#' @param key_match Character. Regular expression used to match against
#'   `names(asp$by_protocol)`. If multiple keys match, the first is used (with a message
#'   when `verbose = TRUE`).
#' @param bl_window Numeric length-2 vector giving the baseline time window in seconds
#'   (e.g. `c(-5, 0)`). Used only if `percent_change` is missing and/or `obj$blMean` is
#'   unavailable/invalid.
#' @param y_mode Character. Y-axis transform for the trace panel:
#'   \describe{
#'     \item{`"percent"`}{plot `% of baseline` (`percent_change`), with reference line at 100}
#'     \item{`"log10_percent"`}{plot `log10(percent_change)` (drops non-positive values), reference at `log10(100)`}
#'     \item{`"signed_log10_delta"`}{plot `sign(%-100) * log10(|%-100| + 1)`, reference at 0}
#'     \item{`"log2_fc"`}{plot `log2((instRate + log_eps_hz)/(baseline + log_eps_hz))`, reference at 0}
#'   }
#' @param y_lim Optional numeric length-2 vector. If provided, sets the visible y-range
#'   via `coord_cartesian(ylim = y_lim)` (x-limits remain the full local time span).
#' @param log_eps_hz Numeric. Small constant added to `instRate` and baseline when
#'   `y_mode = "log2_fc"` to avoid `-Inf` when firing rate approaches 0.
#' @param alpha_trace Numeric in [0, 1]. Alpha for the colored points in the trace panel.
#' @param linewidth Numeric. Line width used for the black trace line in the trace panel.
#' @param heat_dt Numeric. Bin width (seconds) for the heat strip. InstRate is aggregated
#'   within bins and then filled to make a continuous strip.
#' @param ttl_shape Integer. Point shape used for the green TTL down-arrow marker
#'   in the heat panel (default `6`).
#' @param heat_method Character. Aggregation function for instRate within each heat bin:
#'   `"mean"` or `"median"`.
#' @param show_sweep_lines Logical. If `TRUE` and `segment_bounds` is present, draw dotted
#'   vertical lines at sweep/segment boundaries (`segment_bounds$start_time`).
#' @param sweep_line_alpha Numeric in [0, 1]. Alpha for the dotted sweep boundary lines.
#' @param save_png Logical. If `TRUE`, save the combined plot to disk via `ggsave()`.
#' @param png_path Character. Output path for the PNG. If `NULL` and `save_png = TRUE`,
#'   defaults to `file.path(getwd(), paste0(asp$md$cell_name, "-asp-spike-protocol.png"))`.
#' @param fig_w Numeric. Width (in inches) passed to `ggsave()`.
#' @param fig_h Numeric. Height (in inches) passed to `ggsave()`.
#' @param dpi Integer. DPI passed to `ggsave()`.
#' @param verbose Logical. If `TRUE`, print messages about key selection and TTL placement.
#'
#' @return Invisibly returns a list with:
#' \describe{
#'   \item{`protocol_key`}{The protocol key that was plotted.}
#'   \item{`df_trace`}{Filtered/augmented trace data used for plotting (includes `y_plot`).}
#'   \item{`df_heat`}{Binned heat-strip data (time bins and instRate).}
#'   \item{`ttl_times_local`}{The chosen TTL time (closest sampled time to 0).}
#'   \item{`plot`}{The combined patchwork plot object.}
#' }
#'
#' @details
#' If `percent_change` is missing, it is computed as `(instRate / baseline) * 100`,
#' using `obj$blMean` when valid, otherwise recomputed from `bl_window`.
#' For TTL placement in this single-protocol mode, the marker is placed at the time
#' value in `spike_puff_output$time` closest to 0 seconds.
#'
#' @examples
#' \dontrun{
#' # Plot by exact key
#' plot_asp_single_protocol(asp, protocol_key = names(asp$by_protocol)[1])
#'
#' # Plot by regex match
#' plot_asp_single_protocol(asp, key_match = "Standard_Puff", y_mode = "log2_fc")
#'
#' # Set y-limits and save
#' plot_asp_single_protocol(
#'   asp,
#'   key_match = "wash-in",
#'   y_mode = "percent",
#'   y_lim = c(-1, 1),
#'   save_png = TRUE
#' )
#' }
#' @export
plot_asp_single_protocol <- function(asp,
                                     protocol_key = NULL,
                                     key_match = NULL,
                                     bl_window = c(-5, 0),
                                     y_mode = c("percent", "log10_percent", "signed_log10_delta", "log2_fc"),
                                     y_lim = NULL,
                                     x_lim = NULL,
                                     log_eps_hz = 0.05,
                                     alpha_trace = 0.9,
                                     linewidth = 0.35,
                                     heat_dt = 0.10,
                                     ttl_shape = 6,
                                     heat_method = c("mean", "median"),
                                     show_sweep_lines = TRUE,     # uses segment_bounds starts
                                     sweep_line_alpha = 0.35,
                                     save_png = FALSE,
                                     png_path = "figs/exp_single_heat",
                                     save_to_rookery = TRUE,
                                     show_plot = TRUE,
                                     fig_w = 10,
                                     fig_h = 4.5,
                                     dpi = 300,
                                     verbose = TRUE) {

  suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(scales)
    library(tidyr)
    library(tibble)
    library(projHCT)
  })

  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  y_mode <- match.arg(y_mode)
  heat_method <- match.arg(heat_method)

  if (is.null(asp) || is.null(asp$by_protocol) || length(asp$by_protocol) == 0) {
    stop("asp$by_protocol is empty.")
  }
  # -----------------------------
  # helper functions
  # -----------------------------

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

  make_simple_proto_id <- function(cell, proto_index, proto_type, puff, exp) {

    clean <- function(x) {
      x <- as.character(x)
      x[is.na(x) | x == ""] <- "NA"
      x <- trimws(x)
      x <- gsub("[^A-Za-z0-9._-]+", "_", x)
      x
    }

    paste(
      clean(cell),
      sprintf("p%02d", proto_index),
      clean(proto_type),
      clean(puff),
      clean(exp),
      sep = "__"
    )
  }

  # -----------------------------
  # choose protocol (soft-fail)
  # -----------------------------
  keys <- names(asp$by_protocol)

  if (!is.null(protocol_key)) {
    if (!protocol_key %in% keys) {
      if (isTRUE(verbose)) message("protocol_key not found: ", protocol_key, " (returning NULL)")
      return(NULL)
    }
    key <- protocol_key

  } else if (!is.null(key_match)) {
    hit <- keys[grepl(key_match, keys)]

    if (length(hit) == 0) {
      if (isTRUE(verbose)) message("No protocol keys matched key_match='", key_match, "' (returning NULL)")
      return(NULL)
    }

    if (length(hit) > 1 && isTRUE(verbose)) {
      message("Multiple keys matched key_match='", key_match, "'; using first:\n  ", hit[1])
    }
    key <- hit[1]

  } else {
    # default to first (but print options)
    if (isTRUE(verbose)) {
      message("No protocol_key/key_match provided. Available keys:\n- ",
              paste(keys, collapse = "\n- "))
      message("Using first key: ", keys[1])
    }
    key <- keys[1]
  }

  proto_index <- match(key, names(asp$by_protocol))


  proto_type <- dplyr::case_when(
    grepl("^spikeTTL", key) ~ "spikeTTL",
    grepl("^TTLstack", key) ~ "TTLstack",
    TRUE ~ "other"
  )

  obj <- asp$by_protocol[[key]]
  d <- obj$spike_puff_output
  puff <- obj$pro_row$Puff_drug
  exp  <- obj$pro_row$Exp

  file_id <- make_simple_proto_id(
    cell        = asp$md$cell_name,
    proto_index = proto_index,
    proto_type  = proto_type,
    puff        = puff,
    exp         = exp
  )

  if (is.null(d) || !"time" %in% names(d)) stop("spike_puff_output missing or lacks time.")
  if (!"instRate" %in% names(d)) stop("spike_puff_output lacks instRate.")

  # compute percent_change if missing
  if (!"percent_change" %in% names(d)) {
    if (verbose) message("Protocol missing percent_change; computing from instRate.")

    blMean <- obj$blMean
    if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
      bl_keep <- d$time > bl_window[1] & d$time <= bl_window[2]
      blMean <- mean(d$instRate[bl_keep], na.rm = TRUE)
    }
    if (!is.finite(blMean) || blMean <= 0) stop("Cannot compute percent_change: bad baseline mean.")
    d$percent_change <- (d$instRate / blMean) * 100
  }

  d <- d %>%
    dplyr::filter(is.finite(.data$time), is.finite(.data$instRate), is.finite(.data$percent_change)) %>%
    dplyr::arrange(.data$time)

  if (nrow(d) == 0) stop("No usable rows after filtering.")

  # baseline mean for log2_fc (robust)
  blMean <- obj$blMean
  if (is.null(blMean) || !is.finite(blMean) || blMean <= 0) {
    bl_keep <- d$time > bl_window[1] & d$time <= bl_window[2]
    blMean <- mean(d$instRate[bl_keep], na.rm = TRUE)
  }

  # -----------------------------
  # y transform
  # -----------------------------
  if (y_mode == "percent") {
    d$y_plot <- d$percent_change
    y_lab <- "% of baseline instRate"
    y_hline <- 100
  } else if (y_mode == "log10_percent") {
    d <- d %>% dplyr::filter(.data$percent_change > 0)
    d$y_plot <- log10(d$percent_change)
    y_lab <- "log10(% of baseline)"
    y_hline <- log10(100)
  } else if (y_mode == "signed_log10_delta") {
    delta <- d$percent_change - 100
    d$y_plot <- sign(delta) * log10(abs(delta) + 1)
    y_lab <- "sign(%-100) * log10(|%-100|+1)"
    y_hline <- 0
  } else if (y_mode == "log2_fc") {
    # avoid -Inf when instRate goes to 0
    if (is.finite(blMean) && blMean > 0) {
      d$y_plot <- log2((d$instRate + log_eps_hz) / (blMean + log_eps_hz))
    } else {
      # fallback to percent_change (still avoid zeros)
      d <- d %>% dplyr::filter(.data$percent_change > 0)
      d$y_plot <- log2(d$percent_change / 100)
    }
    y_lab <- "log2(instRate / baseline)"
    y_hline <- 0
  }

  if (nrow(d) == 0) stop("No usable rows after y transform.")

  # axis limits
  xlim_local <- range(d$time, na.rm = TRUE)
  rate_lim <- range(d$instRate, na.rm = TRUE)
  max_rate <- max(d$instRate, na.rm = TRUE)

  if (!is.null(y_lim)) {
    ylim_use <- y_lim
  } else {
    yr <- range(d$y_plot, na.rm = TRUE)
    ylim_use <- c(max(yr[1], -1.5), yr[2])
  }


  # -----------------------------
  # TTL time for SINGLE PROTOCOL
  # -----------------------------
  # In this plotting mode, `d$time` is already local with TTL onset at ~0.
  # So do NOT use segment_bounds$start_time (that's concatenated/global).
  ttl_source <- "single_protocol_t0"

  # Use the closest available sample to 0 so the arrow lands exactly on your data grid
  ttl0 <- d$time[which.min(abs(d$time - 0))]
  ttl_times_local <- ttl0

  if (verbose) {
    message("Using protocol key:\n  ", key)
    message("TTL source=", ttl_source, " | ttl_times_local=", round(ttl_times_local, 6))
  }

  df_ttl <- data.frame(
    time = ttl_times_local,
    stringsAsFactors = FALSE
  )

  # sweep boundaries for dotted lines (optional)
  df_sweeps <- NULL
  if (isTRUE(show_sweep_lines) &&
      !is.null(obj$segment_bounds) &&
      all(c("start_time","chunk") %in% names(obj$segment_bounds))) {

    df_sweeps <- obj$segment_bounds %>%
      dplyr::filter(is.finite(.data$start_time)) %>%
      dplyr::distinct(.data$start_time) %>%
      dplyr::arrange(.data$start_time) %>%
      dplyr::rename(x0 = .data$start_time)
  }

  # -----------------------------
  # Heat strip bins (continuous within local span)
  # -----------------------------
  df_heat_raw <- d %>%
    dplyr::mutate(time_bin_start = floor(.data$time / heat_dt) * heat_dt) %>%
    dplyr::group_by(.data$time_bin_start) %>%
    dplyr::summarise(
      instRate = if (heat_method == "mean") mean(.data$instRate, na.rm = TRUE) else median(.data$instRate, na.rm = TRUE),
      .groups = "drop"
    )

  bin_seq <- seq(
    from = floor(xlim_local[1] / heat_dt) * heat_dt,
    to   = floor(xlim_local[2] / heat_dt) * heat_dt,
    by   = heat_dt
  )

  df_heat <- tibble::tibble(time_bin_start = bin_seq) %>%
    dplyr::left_join(df_heat_raw, by = "time_bin_start") %>%
    dplyr::mutate(
      instRate = dplyr::coalesce(.data$instRate, 0),  # <--- key change
      time_bin = .data$time_bin_start + heat_dt/2,
      y = 1
    )


  # -----------------------------
  # Plot: trace (points colored by instRate)
  # -----------------------------

  p_trace_title <- paste(asp$md$cell_name, parse_exp_from_key(key), rank_label(obj), sep = "  |  ")

  p_trace <- ggplot(d, aes(x = time, y = y_plot)) +
    geom_hline(yintercept = y_hline, linewidth = 0.3) +
    geom_line(colour = "black", alpha = 0.25, linewidth = linewidth) +
    geom_point(aes(colour = instRate), alpha = alpha_trace, size = 0.30) +
    scale_colour_gradient(
      low = "navy",          # darker than "blue"
      high = "orange",
      limits = c(0, max_rate),
      oob = scales::squish,
      na.value = "grey60"
    ) +
    guides(colour = "none") +
    labs(
      title = p_trace_title,
      subtitle = key,
      x = NULL,
      y = y_lab
    ) +
    theme_classic() +
    theme(plot.title.position = "plot") +
    scale_x_continuous(limits = xlim_local, expand = c(0, 0)) +
    coord_cartesian(xlim = x_lim %||% xlim_local, ylim = ylim_use, expand = FALSE) +
    theme(plot.margin = margin(t = 5, r = 5, b = 0, l = 5))

  if (!is.null(df_sweeps) && nrow(df_sweeps) > 0) {
    p_trace <- p_trace +
      geom_vline(
        data = df_sweeps,
        aes(xintercept = x0),
        linetype = "dotted",
        linewidth = 0.4,
        alpha = sweep_line_alpha
      )
  }

  # -----------------------------
  # Plot: heat strip + GREEN DOWN ARROWS at TTL
  # -----------------------------
  p_heat <- ggplot(df_heat, aes(x = time_bin, y = y, fill = instRate)) +
    geom_tile(height = 1, width = heat_dt) +
    scale_fill_gradient(
      low = "navy",          # darker than "blue"
      high = "orange",
      limits = c(0, max_rate),
      oob = scales::squish
    ) +
    scale_x_continuous(limits = xlim_local, expand = c(0, 0)) +
    coord_cartesian(xlim = x_lim %||% xlim_local, ylim = c(0.5, 1.8), expand = FALSE) +
    #coord_cartesian(xlim = xlim_local, ylim = c(0.5, 1.8), expand = FALSE) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      plot.margin  = margin(t = 0, r = 5, b = 5, l = 5),
      legend.position = "right"
    ) +
    labs(
      x = "Time (s)",
      y = NULL,
      fill = "instRate (Hz)"
    )

  if (!is.null(df_sweeps) && nrow(df_sweeps) > 0) {
    p_heat <- p_heat +
      geom_vline(
        data = df_sweeps,
        aes(xintercept = x0),
        linetype = "dotted",
        linewidth = 0.4,
        alpha = sweep_line_alpha
      )
  }

  # green arrows (down)
  if (!is.null(df_ttl) && nrow(df_ttl) > 0) {
    p_heat <- p_heat +
      geom_point(
        data = df_ttl,
        aes(x = time, y = 1.62),
        inherit.aes = FALSE,
        shape = ttl_shape,
        size = 2,
        colour = "green4"
      )
  }

  # -----------------------------
  # Combine
  # -----------------------------
  p <- p_trace / p_heat + patchwork::plot_layout(heights = c(3, 1))

  if(isTRUE(show_plot)){
    print(p)
  }




  if (isTRUE(save_to_rookery)) {

    out_dir <- file.path("rookery", asp$md$cell_name)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }

    rookery_png <- file.path(out_dir, paste0(file_id, ".png"))
    rookery_pdf <- file.path(out_dir, paste0(file_id, ".pdf"))

    ggsave(rookery_png, plot = p, width = fig_w, height = fig_h, dpi = dpi)
    ggsave(rookery_pdf, plot = p, width = fig_w, height = fig_h)

    if (isTRUE(verbose)) {
      message("Saved:\n  ", rookery_png, "\n  ", rookery_pdf)
    }
  }

  if (isTRUE(save_png)) {

    is_file <- grepl("\\.png$", png_path, ignore.case = TRUE)

    if (is_file) {
      out_dir <- dirname(png_path)
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
      png_file <- png_path
    } else {
      if (!dir.exists(png_path)) {
        dir.create(png_path, recursive = TRUE, showWarnings = FALSE)
      }
      png_file <- file.path(png_path, paste0(file_id, ".png"))
    }

    ggsave(
      filename = png_file,
      plot = p,
      width = fig_w,
      height = fig_h,
      dpi = dpi
    )

    if (isTRUE(verbose)) {
      message("Saved: ", png_file)
    }
  }


  invisible(list(
    protocol_key = key,
    file_id = file_id,
    df_trace = d,
    df_heat = df_heat,
    ttl_times_local = ttl_times_local,
    plot = p
  ))
}
