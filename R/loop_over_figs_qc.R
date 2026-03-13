make_figs_perkey <- function(
    srt,
    settings,
    heat_dt = 1,
    ttl_shape = 6,
    show_plot = TRUE,
    save_png = TRUE
) {

  y_mode <- settings$y_mode %||% "log2_fc"

  # --- singles: each key can have its own y/x lim ---
  keys <- srt$dfs$asp$protocol_keys

  srt$figs$singles <- setNames(
    lapply(keys, function(key) {

      cfg <- settings$singles[[key]]
      if (is.null(cfg)) cfg <- list()

      plot_asp_single_protocol(
        srt$dfs$asp,
        protocol_key = key,
        heat_dt = heat_dt,
        ttl_shape = ttl_shape,
        y_mode = y_mode,
        y_lim = cfg$y_lim %||% c(NA_real_, NA_real_),
        x_lim = cfg$x_lim %||% NULL,   # only if your function supports it
        show_plot = show_plot,
        save_png = save_png
      )
    }),
    keys
  )

  # --- concat: one shared config ---
  if (length(names(srt$dfs$asp$by_protocol)) >= 2) {
    ccfg <- settings$concat %||% list()
    srt$figs$concat <- plot_asp_concatenated(
      srt$dfs$asp,
      gap_s = 2,
      heat_dt = heat_dt,
      order_by = "wrapped_sweeps",
      y_mode = y_mode,
      y_lim = ccfg$y_lim %||% c(NA_real_, NA_real_),
      x_lim = ccfg$x_lim %||% NULL,   # only if your function supports it
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

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

run_asp_fig_qc_loop_perkey <- function(
    srt,
    settings = NULL,
    show_title = TRUE,
    verbose = TRUE
) {

  parse_lim <- function(txt) {
    txt <- trimws(txt)
    if (!nzchar(txt)) return(NULL)                 # keep
    if (tolower(txt) %in% c("na", "none", "null")) return(c(NA_real_, NA_real_))
    txt2 <- gsub("[\\[\\]()]", "", txt)
    vals <- suppressWarnings(as.numeric(strsplit(txt2, "[, ]+")[[1]]))
    vals <- vals[is.finite(vals)]
    if (length(vals) != 2) return(NA)              # invalid
    vals
  }

  fmt_lim <- function(v) {
    if (is.null(v)) return("NULL")
    paste0("[", v[1], ",", v[2], "]")
  }

  # init settings if missing
  if (is.null(settings)) {
    settings <- list(
      y_mode = "log2_fc",
      singles = list(),      # per-key configs live here
      concat = list(y_lim = c(NA_real_, NA_real_), x_lim = NULL)
    )
  }

  call <- NA_character_

  repeat {

    # (re)make figs using current per-key settings
    srt <- make_figs_perkey(srt, settings = settings)

    title_txt <- if (isTRUE(show_title)) {
      sprintf("%s | y_mode=%s", srt$cell %||% "cell", settings$y_mode %||% "log2_fc")
    } else {
      srt$cell %||% "cell"
    }

    call <- select.list(
      choices = c("keep", "omit", "rerun_figures", "adjust_axes"),
      title   = title_txt
    )

    if (is.null(call) || call == "") next
    if (call %in% c("keep", "omit")) break

    if (call == "adjust_axes") {

      # choose target plot to adjust
      keys <- srt$dfs$asp$protocol_keys
      targets <- c(paste0("single: ", keys), "concat")

      target <- select.list(targets, title = "Which figure do you want to adjust?")
      if (is.null(target) || target == "") next

      # optional: change shared y_mode
      ans <- readline(sprintf("y_mode (shared; current: %s) [Enter keep]: ", settings$y_mode))
      if (nzchar(ans)) settings$y_mode <- ans

      if (startsWith(target, "single: ")) {
        key <- sub("^single: ", "", target)

        cfg <- settings$singles[[key]]
        if (is.null(cfg)) cfg <- list(y_lim = c(NA_real_, NA_real_), x_lim = NULL)

        ans <- readline(sprintf("y_lim for %s 'min,max' or 'na' (current: %s): ",
                                key, fmt_lim(cfg$y_lim)))
        v <- parse_lim(ans)
        if (is.atomic(v) && length(v) == 1 && is.na(v)) {
          if (isTRUE(verbose)) message("Invalid y_lim; keeping previous.")
        } else if (!is.null(v)) {
          cfg$y_lim <- v
        }

        # x_lim optional
        ans <- readline(sprintf("x_lim for %s 'min,max' | 'none' to NULL | Enter keep (current: %s): ",
                                key, fmt_lim(cfg$x_lim)))
        if (nzchar(trimws(ans))) {
          if (tolower(trimws(ans)) %in% c("none","null")) {
            cfg$x_lim <- NULL
          } else {
            v <- parse_lim(ans)
            if (is.atomic(v) && length(v) == 1 && is.na(v)) {
              if (isTRUE(verbose)) message("Invalid x_lim; keeping previous.")
            } else if (!is.null(v)) {
              cfg$x_lim <- v
            }
          }
        }

        settings$singles[[key]] <- cfg

      } else if (target == "concat") {

        ccfg <- settings$concat
        if (is.null(ccfg$y_lim)) ccfg$y_lim <- c(NA_real_, NA_real_)

        ans <- readline(sprintf("y_lim for concat 'min,max' or 'na' (current: %s): ",
                                fmt_lim(ccfg$y_lim)))
        v <- parse_lim(ans)
        if (is.atomic(v) && length(v) == 1 && is.na(v)) {
          if (isTRUE(verbose)) message("Invalid y_lim; keeping previous.")
        } else if (!is.null(v)) {
          ccfg$y_lim <- v
        }

        ans <- readline(sprintf("x_lim for concat 'min,max' | 'none' to NULL | Enter keep (current: %s): ",
                                fmt_lim(ccfg$x_lim)))
        if (nzchar(trimws(ans))) {
          if (tolower(trimws(ans)) %in% c("none","null")) {
            ccfg$x_lim <- NULL
          } else {
            v <- parse_lim(ans)
            if (is.atomic(v) && length(v) == 1 && is.na(v)) {
              if (isTRUE(verbose)) message("Invalid x_lim; keeping previous.")
            } else if (!is.null(v)) {
              ccfg$x_lim <- v
            }
          }
        }

        settings$concat <- ccfg
      }

      next
    }

    # otherwise rerun_figures -> loop continues
  }

  list(
    srt = srt,
    call = call,
    settings = settings
  )
}

