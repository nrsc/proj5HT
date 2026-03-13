# ============================================================
# 5HT child nnest built off HCT header truth
# ============================================================

.local_cache <- new.env(parent = emptyenv())

`%||%` <- function(a, b) if (!is.null(a)) a else b

# What HCT fields are considered "header truth"?
hct_header_fields <- function() c("cell","wd","rd","md","exp","version","files")

hct_header <- function(hct, fields = hct_header_fields()) {
  stopifnot(is.list(hct))
  keep <- intersect(fields, names(hct))
  hct[keep]
}

apply_hct_header <- function(child, header,
                             overwrite = hct_header_fields(),
                             child_files_key = "child_files",
                             preserve_child_files = TRUE) {
  stopifnot(is.list(child), is.list(header))

  saved_child_files <- NULL
  if (preserve_child_files && !is.null(child[[child_files_key]])) {
    saved_child_files <- child[[child_files_key]]
  }

  nm_over <- intersect(overwrite, intersect(names(header), names(child)))
  if (length(nm_over)) child[nm_over] <- header[nm_over]

  nm_add <- setdiff(names(header), names(child))
  if (length(nm_add)) child[nm_add] <- header[nm_add]

  if (preserve_child_files && !is.null(saved_child_files)) {
    child[[child_files_key]] <- saved_child_files
  }

  child
}

parent_state_from_file <- function(path) {
  if (!is.character(path) || length(path) != 1 || !nzchar(path) || !file.exists(path)) return(NULL)
  info <- file.info(path)
  list(
    path  = normalizePath(path, winslash = "/", mustWork = FALSE),
    size  = as.numeric(info$size),
    mtime = as.numeric(info$mtime)
  )
}

parent_state_equal <- function(a, b) {
  if (is.null(a) || is.null(b)) return(FALSE)
  identical(a$path, b$path) &&
    identical(a$size,  b$size) &&
    identical(a$mtime, b$mtime)
}

# Use the same logic as loadHCT "FAST PATH" but return *path* (no readRDS)
find_hct_path_fast <- function(cell) {
  map <- get_idx_map("-hct.rds")
  if (!is.null(map) && cell %in% names(map)) {
    p <- map[[cell]]
    if (file.exists(p)) return(p)
  }
  NA_character_
}



resolve_hct_path <- function(child, cell, source = projHCT::params$rookery) {
  # 1) trust explicit pointer in child (best)
  p <- child$parent$path %||% child$parent$nnest
  if (is.character(p) && length(p) == 1 && file.exists(p)) return(p)

  # 2) fast index lookup (no file read)
  p2 <- find_hct_path_fast(cell)
  if (!is.na(p2) && file.exists(p2)) return(p2)

  # 3) optional fallback disk search (slower, but still no readRDS)
  cell_pat <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", cell)
  pat <- paste0(cell_pat, "-hct\\.rds$")
  lf <- list.files(source, pattern = pat, recursive = TRUE, full.names = TRUE)
  lf <- unique(lf[file.exists(lf)])
  if (!length(lf)) return(NA_character_)
  lf[1]
}


get_parent_path_from_child <- function(child, cell) {
  # prefer the authoritative pointer you already store
  p <- child$parent$path %||% child$parent$nnest
  if (!is.null(p) && is.character(p) && length(p) == 1) return(p)

  # fallback: infer from rd/cell
  if (!is.null(child$rd) && !is.null(cell)) {
    return(file.path(child$rd, paste0(cell, "-hct.rds")))
  }
  NA_character_
}

# returns the on-disk path for a given cell+tag without readRDS()
find_child_path_fast <- function(cell, tag, source = projHCT::params$rookery) {

  tag <- normalize_tag(tag)

  # ---- FAST PATH: cached index lookup ----
  map <- get_idx_map(tag)
  if (!is.null(map) && cell %in% names(map)) {
    p <- map[[cell]]
    if (file.exists(p)) return(p)
  }

  # ---- FALLBACK: disk scan (slow, rarely used) ----
  if (!is.character(source) || !length(source) || !dir.exists(source)) {
    return(NA_character_)
  }

  cell_pat <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", cell)
  pat <- paste0(cell_pat, tag, "$")

  lf <- list.files(source, pattern = pat, recursive = TRUE, full.names = TRUE)
  lf <- unique(lf[file.exists(lf)])
  if (!length(lf)) return(NA_character_)

  lf[1]
}

## Normalize the tag
normalize_tag <- function(tag) {
  if (!grepl("\\.rds$", tag, ignore.case = TRUE)) tag <- paste0(tag, ".rds")
  if (!startsWith(tag, "-")) tag <- paste0("-", tag)
  tag
}

get_idx_map <- function(tag) {
  key <- paste0("idx_", tag)
  if (exists(key, envir = .local_cache, inherits = FALSE)) {
    return(get(key, envir = .local_cache))
  }

  idx <- switch(
    tag,
    "-srt.rds" = projHCT::sheets$files$nnest_srt_files,
    "-osc.rds" = projHCT::sheets$files$nnest_osc_files,
    "-hct.rds" = projHCT::sheets$files$nnest_files,
    NULL
  )

  if (is.null(idx)) return(NULL)

  map <- setNames(idx, vapply(idx, nphys::fileD, character(1)))
  assign(key, map, envir = .local_cache)
  map
}


#' Create / update the proj5HT child payload nnest for a cell (srt)
#'
#' HCT (-hct.rds) is the authoritative header (static state).
#' proj5HT child (-srt.rds) is mutable analysis payload (dfs/ana/etc.).
#'
#' @param x cell name (character) OR an HCT nnest/list with $cell
#' @param tag_child filename tag for child payload (default "-srt.rds")
#' @param save_nnest write payload to disk
#' @param update if TRUE, preserves existing payload (dfs/ana/etc.) when rebuilding
#' @param write if FALSE, dry run (no file writes)
#' @param verbose messages
#' @param payload_keys which top-level keys to carry forward from an existing payload on update
#' @param mimic_header if TRUE, attach header-like fields (cell/wd/rd/md/files) onto returned child
#'
#' @return child nnest object (list). Status info is attached as attr(child, "status")
#' @export
nnest5HT <- function(
    x,
    tag_child = "-srt.rds",
    save_nnest = TRUE,
    update = TRUE,
    write = TRUE,
    verbose = TRUE,
    payload_keys = c("dfs", "ana", "smry", "params", "results", "log"),
    mimic_header = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(...)

  # ---- resolve cell name ----
  cell <- if (is.list(x)) x$cell else x
  if (!is.character(cell) || length(cell) != 1 || !nzchar(cell)) {
    stop("nnest5HT: x must be a single cell name or an nnest/list with $cell")
  }

  # ---- normalize tag_child safely ----
  tag_child <- normalize_tag(tag_child)

  # ---- load / ensure parent HCT header truth ----
  hct <- loadHCT(cell, tag = "-hct.rds")
  if (is.null(hct)) {
    msg("No HCT header found for ", cell, " — building via nnestHCT()")
    hct <- nnestHCT(cell)
    if (is.null(hct)) stop("failed_to_build_hct for cell: ", cell)
  }

  # header (static truth)
  header <- hct_header(hct)
  header$files <- header$files %||% list()

  # canonical header nnest path
  header_nnest_path <- file.path(header$rd, paste0(cell, "-hct.rds"))
  header$files$nnest <- header$files$nnest %||% header_nnest_path

  # ---- define where the child payload lives ----
  child_path <- file.path(header$rd, paste0(cell, tag_child))

  # ---- carry forward payload from existing child ----
  carry <- list()
  if (isTRUE(update) && file.exists(child_path)) {
    old <- tryCatch(readRDS(child_path), error = function(e) NULL)
    if (is.list(old)) {
      for (k in intersect(payload_keys, names(old))) carry[[k]] <- old[[k]]
    }
  }

  # ---- build the child payload object ----
  parent_path <- header$files$nnest %||% NA_character_

  child <- list(
    project = "5HT",
    module  = "srt",

    parent  = list(
      cell = cell,
      tag  = "-hct.rds",
      path = parent_path
    ),

    parent_state = if (!is.na(parent_path) && nzchar(parent_path) && file.exists(parent_path)) {
      parent_state_from_file(parent_path)
    } else {
      NULL
    },

    child_files = list(
      nnest = child_path,
      md_snapshot = file.path(header$rd, paste0(cell, "-5ht_md_snapshot.csv"))
      # smry        = file.path(header$rd, paste0(cell, "-5ht_smry.csv"))
    )
  )

  if (length(carry)) child <- c(child, carry)

  # ---- optionally mimic header fields for downstream convenience ----
  if (isTRUE(mimic_header)) {
    child$cell  <- header$cell
    child$wd    <- header$wd
    child$rd    <- header$rd
    child$md    <- header$md
    child$files <- header$files
  }

  # ---- write snapshot + payload (with guards) ----
  if (isTRUE(write) && isTRUE(save_nnest)) {
    if (!dir.exists(header$rd)) dir.create(header$rd, recursive = TRUE, showWarnings = FALSE)

    # Guard md_snapshot path
    md_out <- child$child_files$md_snapshot
    if (length(md_out) == 1 && nzchar(md_out)) {
      tryCatch(
        utils::write.csv(header$md, file = md_out, row.names = FALSE),
        error = function(e) msg("Could not write md snapshot for ", cell, ": ", conditionMessage(e))
      )
    } else {
      msg("Skipping md snapshot (invalid path) for ", cell)
    }

    # Guard child nnest path
    out_path <- child$child_files$nnest
    if (length(out_path) != 1 || !nzchar(out_path)) {
      stop("Invalid child_files$nnest path for cell ", cell)
    }

    saveRDS(child, file = out_path)
    msg("Saved 5HT payload: ", out_path)

  } else if (isTRUE(save_nnest)) {
    msg("[DRY RUN] Would save 5HT payload: ", child$child_files$nnest)
  }

  # attach status info as attribute (optional)
  attr(child, "status") <- list(cell = cell, ok = TRUE, payload_path = child$child_files$nnest)

  child
}

#' Load proj5HT child nnest (e.g. srt) and rehydrate header truth from HCT
#'
#' @param x cell name OR an nnest-like object OR a status list with $payload_path
#' @param tag child tag, e.g. "-srt.rds"
#' @param rehydrate logical; if TRUE, inject HCT header truth into child in-memory
#' @param overwrite header fields to overwrite in child when rehydrating
#' @param verbose logical
#'
#' @return child nnest object (full structure), or NULL if not found
#' @export
load5HT <- function(
    x,
    tag = "-srt.rds",
    rehydrate = TRUE,
    overwrite = c("cell","wd","rd","md","exp","version","files"),
    verbose = TRUE
) {
  msg <- function(...) if (isTRUE(verbose)) message(...)

  # normalize tag
  tag <- normalize_tag(tag)

  # ---- helper: infer cell name from inputs ----
  infer_cell <- function(x) {
    if (is.character(x) && length(x) == 1) return(x)
    if (is.list(x) && !is.null(x$cell) && is.character(x$cell) && length(x$cell) == 1) return(x$cell)
    if (is.list(x) && !is.null(x$payload_path) && is.character(x$payload_path) && length(x$payload_path) == 1) {
      # try to parse cell from file name if you have nphys::fileD
      return(tryCatch(nphys::fileD(x$payload_path), error = function(e) NA_character_))
    }
    NA_character_
  }

  #cell <- infer_cell(x)
  cell <- infer_cell(cell)

  # ---- FAST PATH: cell name + no rehydrate ----
  if (isFALSE(rehydrate) &&
      is.character(cell) &&
      length(cell) == 1 &&
      nzchar(cell)) {

    child_path <- find_child_path_fast(cell, tag)

    if (!is.character(child_path) || length(child_path) != 1 || !file.exists(child_path)) {
      if (isTRUE(verbose)) message("No child file found: ", child_path)
      child = nnest5HT(cell)
      return(child)
    }

    return(readRDS(child_path))

  }

  # NOTE:
  # The code below supports legacy / flexible inputs:
  # - nnest-like lists
  # - status objects with $payload_path
  # - rehydration against HCT header truth
  #
  # In practice, proj5HT workflows overwhelmingly call:
  #   load5HT(cell, tag = "-srt.rds", rehydrate = FALSE)
  #
  # Do not optimize the paths below at the expense of the fast path above.


  # ---- if x is a status list with payload_path, load that directly ----
  if (is.list(x) && !is.null(x$payload_path) && is.character(x$payload_path) && length(x$payload_path) == 1) {
    if (!file.exists(x$payload_path)) {
      msg("payload_path does not exist: ", x$payload_path)
      return(NULL)
    }
    child <- readRDS(x$payload_path)
    # if we couldn't infer cell, try from child
    if (is.na(cell) || !nzchar(cell)) cell <- infer_cell(child)
  } else {
    # ---- normal path: load header to determine rd, then load child by tag ----
    if (is.na(cell) || !nzchar(cell)) {
      stop("load5HT: could not determine cell name from x")
    }

    # hct <- loadHCT(cell, tag = "-hct.rds")
    # if (is.null(hct)) {
    #   msg("No HCT header found for ", cell)
    #   return(NULL)
    # }

    child_path <- find_child_path_fast(cell, tag)
    if (!file.exists(child_path)) {
      msg("No child file found: ", child_path)
      return(NULL)
    }
    child <- readRDS(child_path)
  }

  # ---- rehydrate header truth into child (recommended) ----
  if (isTRUE(rehydrate)) {

    # --- fast gate: check parent hct file state WITHOUT loading the header object ---
    parent_path <- resolve_hct_path(child, cell)

    if (!is.na(parent_path) && file.exists(parent_path)) {
      current_state <- parent_state_from_file(parent_path)
      stored_state  <- child$parent_state %||% NULL

      if (parent_state_equal(stored_state, current_state)) {
        # Parent header truth unchanged → return child as-is (no loadHCT, no overwrite)
        return(child)
      }
    }

    # --- slow path: parent changed or unknown → load & rehydrate ---
    hct <- loadHCT(cell, tag = "-hct.rds")
    if (!is.null(hct)) {

      nm <- intersect(overwrite, intersect(names(child), names(hct)))
      if (length(nm)) child[nm] <- hct[nm]

      missing <- setdiff(intersect(overwrite, names(hct)), names(child))
      if (length(missing)) child[missing] <- hct[missing]

      child$child_files <- child$child_files %||% list()
      child$child_files$nnest <- file.path(hct$rd, paste0(cell, tag))

      # ---- parent pointers: reuse parent_path; fall back to canonical hct path ----
      canonical_parent_path <- file.path(hct$rd, paste0(cell, "-hct.rds"))
      if (is.na(parent_path) || !nzchar(parent_path) || !file.exists(parent_path)) {
        parent_path <- canonical_parent_path
      }

      child$parent <- child$parent %||% list()
      child$parent$path  <- parent_path
      child$parent$nnest <- parent_path
      child$parent$tag   <- "-hct.rds"

      child$project <- child$project %||% "5HT"
      child$module  <- child$module  %||% "srt"

      # refresh fingerprint using the same parent_path (no extra resolve_hct_path call)
      if (!is.na(parent_path) && file.exists(parent_path)) {
        child$parent_state <- parent_state_from_file(parent_path)
      }
    }

    # IMPORTANT: ensure this block evaluates to the nnest, not parent_state
    child
  }

  child

}

