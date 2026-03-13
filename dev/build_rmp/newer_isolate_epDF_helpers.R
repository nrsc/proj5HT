####

#' Is an epochDF usable after cleaning?
#'
#' @param epDF Raw epochDF.
#' @param t0_col Start time column name.
#' @param t1_col End time column name.
#'
#' @return Logical scalar.
#' @keywords internal
is_epochDF_usable <- function(epDF, t0_col = "T0", t1_col = "T1") {
  ep <- clean_epochDF(epDF, t0_col = t0_col, t1_col = t1_col)
  !is.null(ep) && nrow(ep) > 0
}
#' Extract stimulus vector for a single sweep from an NWB file
#'
#' Lightweight wrapper around \code{nphys::extractNWB()} that safely
#' retrieves the stimulus trace for a given sweep.
#'
#' @param nwb An NWB object returned by \code{nphys::extractNWB()}.
#' @param nm Sweep name (character), e.g. "sweep_060".
#'
#' @return Numeric vector of stimulus values, or \code{NULL} if unavailable.
#'
#' @keywords internal
get_stim_for_sweep <- function(nwb, nm) {

  # Defensive checks
  if (is.null(nwb) || is.null(nm)) return(NULL)
  if (is.null(nwb$nwb_file)) return(NULL)

  stim <- tryCatch(
    nphys::extractNWB(
      nwb$nwb_file,
      stimulus_sweeps = nm,
      slim = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(stim)) return(NULL)

  as.numeric(stim)
}
