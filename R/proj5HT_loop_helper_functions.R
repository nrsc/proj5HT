#### proj5HT Helper functions
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
# helper to (re)make figures
make_figs <- function(srt, y_lim_single = NULL, y_lim_concat = NULL, x_lim_by_protocol = c(-10, 50), y_mode) {

  srt$figs$singles <- setNames(
    lapply(srt$dfs$asp$protocol_keys, function(key) {
      plot_asp_single_protocol(
        srt$dfs$asp,
        protocol_key = key,
        heat_dt = 1,
        ttl_shape = 6,
        y_mode = y_mode,
        y_lim = y_lim_single,
        x_lim = x_lim_by_protocol,
        show_plot = TRUE,
        save_png = TRUE
      )
    }),
    srt$dfs$asp$protocol_keys
  )

  if (length(names(srt$dfs$asp$by_protocol)) >= 2) {
    srt$figs$concat <- plot_asp_concatenated(
      srt$dfs$asp,
      gap_s = 2,
      heat_dt = 1,
      order_by = "wrapped_sweeps",
      y_mode = y_mode,
      y_lim = y_lim_concat,
      x_lim_by_protocol = x_lim_by_protocol,
      show_gap_lines = TRUE,
      save_png = TRUE,
      png_path = "figs/exp_concat_heat/",
      verbose = TRUE
    )
  } else {
    srt$figs$concat <- NULL
  }

  srt$figs$washin <- plot_sidebyside_heat_wash(srt, show_plot = TRUE)

  srt
}


#' Interactive QC loop for ASP figures with adjustable axis settings
#'
#' Repeatedly generates per-protocol and concatenated ASP figures, displays them,
#' and prompts the user to choose keep/omit/rerun, with an option to adjust axis
#' settings (y-mode, y-lims, optional x-lims) between reruns.
#'
#' @param srt A nnest-like object containing `srt$dfs$asp` and `srt$cell`.
#' @param make_figs_fn Function that takes `(srt, y_lim_single, y_lim_concat, y_mode, x_lim_single, x_lim_concat)`
#'   and returns an updated `srt` with figures stored (e.g. `srt$figs$singles`, `srt$figs$concat`, `srt$figs$washin`).
#' @param y_mode Character scalar. Initial y-mode to pass through to plotting.
#' @param y_lim_single Numeric length-2. Initial y-lims for single-protocol plots.
#' @param y_lim_concat Numeric length-2. Initial y-lims for concatenated plots.
#' @param x_lim_single NULL or numeric length-2. Optional x-lims for single-protocol plots.
#' @param x_lim_concat NULL or numeric length-2. Optional x-lims for concatenated plots.
#' @param show_title Logical. If TRUE, show current settings in the select.list title.
#' @param verbose Logical. If TRUE, prints small messages (e.g. invalid inputs).
#'
#' @return A list with:
#' \describe{
#'   \item{`srt`}{Updated `srt` with figures attached.}
#'   \item{`call`}{Final decision: "keep" or "omit".}
#'   \item{`settings`}{List of final settings used (y_mode, y_lim_single, y_lim_concat, x_lim_single, x_lim_concat).}
#' }
#' @export
run_asp_fig_qc_loop <- function(
    srt,
    make_figs_fn,
    y_mode = "log2_fc",
    y_lim_single = c(NA_real_, NA_real_),
    y_lim_concat = c(NA_real_, NA_real_),
    x_lim_single = NULL,
    x_lim_concat = NULL,
    show_title = TRUE,
    verbose = TRUE
) {

  stopifnot(is.function(make_figs_fn))

  parse_lim <- function(txt) {
    txt <- trimws(txt)
    if (!nzchar(txt)) return(NULL) # means "no change"
    if (tolower(txt) %in% c("na", "none", "null")) return(c(NA_real_, NA_real_))
    txt2 <- gsub("[\\[\\]()]", "", txt)
    vals <- suppressWarnings(as.numeric(strsplit(txt2, "[, ]+")[[1]]))
    vals <- vals[is.finite(vals)]
    if (length(vals) != 2) return(NA) # invalid
    vals
  }

  fmt_lim <- function(v) {
    if (is.null(v)) return("NULL")
    if (length(v) != 2) return("<?>")
    paste0("[", v[1], ",", v[2], "]")
  }

  call <- NA_character_

  repeat {

    # (re)make figs with current settings
    srt <- make_figs_fn(
      srt,
      y_lim_single = y_lim_single,
      y_lim_concat = y_lim_concat,
      y_mode       = y_mode,
      x_lim_single = x_lim_single,
      x_lim_concat = x_lim_concat
    )

    title_txt <- if (isTRUE(show_title)) {
      sprintf(
        "%s | y_mode=%s | y_single=%s | y_concat=%s | x_single=%s | x_concat=%s",
        srt$cell %||% "cell",
        y_mode,
        fmt_lim(y_lim_single),
        fmt_lim(y_lim_concat),
        fmt_lim(x_lim_single),
        fmt_lim(x_lim_concat)
      )
    } else {
      srt$cell %||% "cell"
    }

    call <- select.list(
      choices = c("keep", "omit", "rerun_figures", "adjust_axes"),
      title   = title_txt
    )

    # user canceled -> ask again
    if (is.null(call) || identical(call, "")) next

    if (call %in% c("keep", "omit")) break

    if (call == "adjust_axes") {

      cat("\n--- Adjust axes (press Enter to keep current) ---\n")

      ans <- readline(sprintf("y_mode (current: %s): ", y_mode))
      if (nzchar(ans)) y_mode <- ans

      ans <- readline(sprintf("y_lim_single 'min,max' or 'na' (current: %s): ", fmt_lim(y_lim_single)))
      v <- parse_lim(ans)
      if (is.atomic(v) && length(v) == 1 && is.na(v)) {
        if (isTRUE(verbose)) message("Invalid y_lim_single; keeping previous.")
      } else if (!is.null(v)) {
        y_lim_single <- v
      }

      ans <- readline(sprintf("y_lim_concat 'min,max' or 'na' (current: %s): ", fmt_lim(y_lim_concat)))
      v <- parse_lim(ans)
      if (is.atomic(v) && length(v) == 1 && is.na(v)) {
        if (isTRUE(verbose)) message("Invalid y_lim_concat; keeping previous.")
      } else if (!is.null(v)) {
        y_lim_concat <- v
      }

      ans <- readline(sprintf("x_lim_single 'min,max' | 'none' to NULL | Enter keep (current: %s): ", fmt_lim(x_lim_single)))
      if (nzchar(trimws(ans))) {
        if (tolower(trimws(ans)) %in% c("none", "null")) {
          x_lim_single <- NULL
        } else {
          v <- parse_lim(ans)
          if (is.atomic(v) && length(v) == 1 && is.na(v)) {
            if (isTRUE(verbose)) message("Invalid x_lim_single; keeping previous.")
          } else if (!is.null(v)) {
            x_lim_single <- v
          }
        }
      }

      ans <- readline(sprintf("x_lim_concat 'min,max' | 'none' to NULL | Enter keep (current: %s): ", fmt_lim(x_lim_concat)))
      if (nzchar(trimws(ans))) {
        if (tolower(trimws(ans)) %in% c("none", "null")) {
          x_lim_concat <- NULL
        } else {
          v <- parse_lim(ans)
          if (is.atomic(v) && length(v) == 1 && is.na(v)) {
            if (isTRUE(verbose)) message("Invalid x_lim_concat; keeping previous.")
          } else if (!is.null(v)) {
            x_lim_concat <- v
          }
        }
      }

      next
    }

    # otherwise: rerun_figures -> loop continues with same settings
  }

  # small infix helper (avoid importing rlang/purrr here)
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

  list(
    srt = srt,
    call = call,
    settings = list(
      y_mode = y_mode,
      y_lim_single = y_lim_single,
      y_lim_concat = y_lim_concat,
      x_lim_single = x_lim_single,
      x_lim_concat = x_lim_concat
    )
  )
}
