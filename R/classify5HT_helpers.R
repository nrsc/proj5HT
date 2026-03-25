# ============================================================
# classify5HT_helpers.R — Internal helpers for 5-HT response classifier
# ============================================================
# These are non-exported utility functions used by classify5HT.R
# and plot_classify5HT.R.
# ============================================================

#' Assert that a data.frame has required columns
#' @param df data.frame
#' @param cols character vector of required column names
#' @param df_name label used in error message
#' @keywords internal
assert_has_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) {
    stop(df_name, " is missing required columns: ", paste(miss, collapse = ", "))
  }

  invisible(TRUE)
}

#' Coalesce character vector: replace NA/empty with default
#' @keywords internal
coalesce_chr <- function(x, default = NA_character_) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- default
  x
}

#' Return the first non-missing entry in a vector
#' @keywords internal
first_non_missing <- function(x, default = NA) {
  idx <- which(!is.na(x) & x != "")
  if (length(idx) == 0) return(default)
  x[idx[1]]
}

#' Smooth a numeric vector via rollmean, rollmedian, or loess
#'
#' @param y numeric vector to smooth
#' @param method one of "rollmean", "rollmedian", "loess"
#' @param k window width for rolling methods
#' @param loess_span span for loess
#' @param x time axis (required for loess)
#' @return numeric vector of smoothed values
#' @keywords internal
smooth_vec <- function(y,
                       method = c("rollmean", "rollmedian", "loess"),
                       k = 3,
                       loess_span = 0.2,
                       x = NULL) {
  method <- match.arg(method)

  if (length(y) == 0) return(y)

  out <- switch(
    method,
    rollmean   = zoo::rollmean(y, k = k, fill = NA, align = "center"),
    rollmedian = zoo::rollmedian(y, k = k, fill = NA, align = "center"),
    loess = {
      if (is.null(x)) stop("smooth_vec(method='loess') requires x")
      ok <- is.finite(x) & is.finite(y)
      if (sum(ok) < 6) {
        rep(NA_real_, length(y))
      } else {
        fit <- stats::loess(y[ok] ~ x[ok], span = loess_span, degree = 1)
        tmp <- rep(NA_real_, length(y))
        tmp[ok] <- stats::predict(fit, newdata = x[ok])
        tmp
      }
    }
  )

  as.numeric(out)
}

#' Find contiguous TRUE runs in a logical vector (with time positions)
#' @param time numeric time axis
#' @param flag logical vector (same length as time)
#' @return tibble with start, end, dur for each TRUE run
#' @keywords internal
find_runs <- function(time, flag) {
  stopifnot(length(time) == length(flag))

  if (length(flag) == 0) {
    return(tibble::tibble(start = numeric(0), end = numeric(0), dur = numeric(0)))
  }

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  idx <- which(r$values)

  if (length(idx) == 0) {
    return(tibble::tibble(start = numeric(0), end = numeric(0), dur = numeric(0)))
  }

  tibble::tibble(
    start = time[starts[idx]],
    end   = time[ends[idx]],
    dur   = time[ends[idx]] - time[starts[idx]]
  )
}

#' Does a persistent excursion exist in a time window?
#'
#' Tests whether the smoothed trace exceeds a threshold continuously for
#' at least \code{min_persist_s} seconds inside the given window.
#'
#' @param time numeric time vector
#' @param y numeric smoothed trace (same length)
#' @param window numeric length-2 (start, end)
#' @param direction "exc" or "inh"
#' @param baseline reference level (default 100)
#' @param thr threshold in percent-of-baseline units
#' @param min_persist_s minimum duration (seconds) for a valid run
#' @return logical
#' @keywords internal
has_persistent_excursion <- function(time,
                                     y,
                                     window,
                                     direction = c("exc", "inh"),
                                     baseline = 100,
                                     thr = 10,
                                     min_persist_s = 0.5) {
  direction <- match.arg(direction)

  keep <- is.finite(time) & is.finite(y) &
    time >= window[1] & time <= window[2]

  time0 <- time[keep]
  y0    <- y[keep]

  if (length(time0) < 2) return(FALSE)

  flag <- if (direction == "exc") {
    y0 >= (baseline + thr)
  } else {
    y0 <= (baseline - thr)
  }

  if (!any(flag)) return(FALSE)

  runs <- find_runs(time0, flag)
  if (nrow(runs) == 0) return(FALSE)

  any(runs$dur >= min_persist_s, na.rm = TRUE)
}

#' Get the extreme (min or max) value inside a time window
#' @keywords internal
get_extreme_in_window <- function(time, y, window, which = c("min", "max")) {
  which <- match.arg(which)

  keep <- is.finite(time) & is.finite(y) &
    time >= window[1] & time <= window[2]

  time0 <- time[keep]
  y0    <- y[keep]

  if (length(y0) == 0) {
    return(list(time = NA_real_, value = NA_real_))
  }

  idx <- if (which == "min") which.min(y0) else which.max(y0)
  list(time = time0[idx], value = y0[idx])
}

#' Bin a magnitude value into mild / moderate / strong / none
#' @keywords internal
classify_magnitude <- function(mag, mild = 10, moderate = 25, strong = 50) {
  dplyr::case_when(
    !is.finite(mag) ~ NA_character_,
    mag >= strong   ~ "strong",
    mag >= moderate ~ "moderate",
    mag >= mild     ~ "mild",
    TRUE            ~ "none"
  )
}

#' Standard colour palette for response classes
#' @keywords internal
response_cols <- c(
  strong_inhibition   = "#B2182B",
  moderate_inhibition = "#D6604D",
  mild_inhibition     = "#F4A582",
  no_change           = "grey80",
  mild_excitation     = "#92C5DE",
  moderate_excitation = "#4393C3",
  strong_excitation   = "#2166AC"
)

#' Pretty-print labels for response classes
#' @keywords internal
pretty_labels <- c(
  strong_inhibition   = "Strong inhibition",
  moderate_inhibition = "Moderate inhibition",
  mild_inhibition     = "Mild inhibition",
  no_change           = "No change",
  mild_excitation     = "Mild excitation",
  moderate_excitation = "Moderate excitation",
  strong_excitation   = "Strong excitation"
)
