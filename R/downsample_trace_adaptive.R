# ------------------------------------------------------------
# Adaptive downsample voltage trace (preserves APs)
# ------------------------------------------------------------

#' Adaptive Downsampling
#'
#' @param v
#' @param sr
#' @param dvdt_thresh_mV_ms
#' @param slow_bin_ms
#' @param fast_keep_every
#' @param max_points
#' @param thin_slow_every
#'
#' @returns
#' @export
#'
#' @examples
downsample_trace_adaptive <- function(v,
                                      sr,
                                      # threshold for "fast" regions, in mV/ms
                                      dvdt_thresh_mV_ms = 5,
                                      # for slow regions: bin size in milliseconds
                                      slow_bin_ms = 10,
                                      # for fast regions: keep every Nth point (1 = keep all)
                                      fast_keep_every = 1,
                                      # hard cap on output points (NULL = no cap)
                                      max_points = 50000,
                                      # if max_points is exceeded: extra thinning for slow-only points
                                      thin_slow_every = 2) {

  stopifnot(is.numeric(v), length(v) >= 5, is.numeric(sr), length(sr) == 1, sr > 0)

  n <- length(v)
  dt <- 1 / sr

  # time vector (seconds)
  t <- (seq_len(n) - 1) * dt

  # derivative in mV/ms
  dv <- c(0, diff(v))
  dvdt_mV_ms <- (dv / dt) / 1000

  fast <- abs(dvdt_mV_ms) >= dvdt_thresh_mV_ms
  fast[1] <- fast[n] <- TRUE  # keep endpoints

  # helper: return indices for min+max within a window [i:j]
  idx_minmax <- function(i, j) {
    seg <- v[i:j]
    kmin <- which.min(seg)
    kmax <- which.max(seg)
    out <- sort(unique(c(i + kmin - 1, i + kmax - 1)))
    out
  }

  # walk through trace and build kept indices
  keep_idx <- integer(0)

  i <- 1L
  while (i <= n) {
    if (fast[i]) {
      # FAST segment: keep densely
      j <- i
      while (j < n && fast[j + 1L]) j <- j + 1L

      seg_idx <- i:j
      if (fast_keep_every > 1L) seg_idx <- seg_idx[seq(1, length(seg_idx), by = fast_keep_every)]

      keep_idx <- c(keep_idx, seg_idx)
      i <- j + 1L
    } else {
      # SLOW segment: compress with min/max per bin
      j <- i
      while (j < n && !fast[j + 1L]) j <- j + 1L

      # binning in samples
      bin_n <- max(2L, as.integer(round((slow_bin_ms / 1000) * sr)))
      starts <- seq(i, j, by = bin_n)

      seg_keep <- integer(0)
      for (s in starts) {
        e <- min(j, s + bin_n - 1L)
        seg_keep <- c(seg_keep, idx_minmax(s, e))
      }

      keep_idx <- c(keep_idx, seg_keep)
      i <- j + 1L
    }
  }

  keep_idx <- sort(unique(keep_idx))

  # enforce endpoints
  keep_idx <- sort(unique(c(1L, keep_idx, n)))

  # If too big: thin only the SLOW points further (do NOT thin fast points)
  if (!is.null(max_points) && length(keep_idx) > max_points) {
    is_fast_kept <- fast[keep_idx]

    fast_idx <- keep_idx[is_fast_kept]
    slow_idx <- keep_idx[!is_fast_kept]

    # thin slow points
    if (length(slow_idx) > 0) {
      slow_idx <- slow_idx[seq(1, length(slow_idx), by = thin_slow_every)]
    }

    keep_idx <- sort(unique(c(fast_idx, slow_idx, 1L, n)))

    # still too big? increase thinning (slow only) until under cap
    while (!is.null(max_points) && length(keep_idx) > max_points && length(slow_idx) > 10) {
      thin_slow_every <- thin_slow_every + 1L
      slow_idx2 <- slow_idx[seq(1, length(slow_idx), by = thin_slow_every)]
      keep_idx <- sort(unique(c(fast_idx, slow_idx2, 1L, n)))
      slow_idx <- slow_idx2
    }
  }

  data.frame(
    idx = keep_idx,
    t = t[keep_idx],
    v = v[keep_idx],
    dvdt_mV_ms = dvdt_mV_ms[keep_idx],
    fast = fast[keep_idx],
    stringsAsFactors = FALSE
  )
}
