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

