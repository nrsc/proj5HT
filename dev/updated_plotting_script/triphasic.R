# Find contiguous TRUE runs for a logical vector and return their time spans
.find_runs <- function(time, flag) {
  stopifnot(length(time) == length(flag))
  if (length(flag) == 0) return(tibble::tibble(start = numeric(0), end = numeric(0), dur = numeric(0)))

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

# Triphasic detector: excitation -> inhibition -> excitation (E->I->E)
detect_triphasic_exc_inh_exc <- function(d,
                                         time_col = "time",
                                         y_col = "pct_smooth",   # use smoothed
                                         window_total = c(0, 20),# scan window
                                         thr = 10,               # min_excursion_pct
                                         min_persist_s = 0.5,    # persistence for each phase
                                         min_gap_s = 0.05,       # avoid immediately-adjacent flips
                                         baseline = 100) {

  stopifnot(is.data.frame(d), time_col %in% names(d), y_col %in% names(d))

  d0 <- d %>%
    dplyr::mutate(
      .time = as.numeric(.data[[time_col]]),
      .y    = as.numeric(.data[[y_col]])
    ) %>%
    dplyr::filter(is.finite(.data$.time), is.finite(.data$.y),
                  .data$.time >= window_total[1], .data$.time <= window_total[2]) %>%
    dplyr::arrange(.data$.time)

  if (nrow(d0) < 5) {
    return(list(
      is_triphasic = FALSE,
      triphasic_type = NA_character_,
      e1_start = NA_real_, e1_end = NA_real_,
      i_start  = NA_real_, i_end  = NA_real_,
      e2_start = NA_real_, e2_end = NA_real_,
      reason = "too_few_points"
    ))
  }

  exc_flag <- d0$.y >= (baseline + thr)
  inh_flag <- d0$.y <= (baseline - thr)

  exc_runs <- .find_runs(d0$.time, exc_flag) %>% dplyr::filter(.data$dur >= min_persist_s)
  inh_runs <- .find_runs(d0$.time, inh_flag) %>% dplyr::filter(.data$dur >= min_persist_s)

  if (nrow(exc_runs) == 0 || nrow(inh_runs) == 0) {
    return(list(
      is_triphasic = FALSE,
      triphasic_type = NA_character_,
      e1_start = NA_real_, e1_end = NA_real_,
      i_start  = NA_real_, i_end  = NA_real_,
      e2_start = NA_real_, e2_end = NA_real_,
      reason = "missing_required_runs"
    ))
  }

  # We want E1 then I then E2 in time order, with small gaps allowed
  # Choose the earliest valid E1, then earliest I after E1, then earliest E2 after I.
  for (i_e1 in seq_len(nrow(exc_runs))) {
    e1 <- exc_runs[i_e1, ]

    # require I starts after E1 ends (+gap)
    inh_after <- inh_runs %>% dplyr::filter(.data$start >= (e1$end + min_gap_s))
    if (nrow(inh_after) == 0) next

    for (i_i in seq_len(nrow(inh_after))) {
      inh <- inh_after[i_i, ]

      # require E2 starts after I ends (+gap)
      exc_after <- exc_runs %>% dplyr::filter(.data$start >= (inh$end + min_gap_s))
      if (nrow(exc_after) == 0) next

      e2 <- exc_after[1, ]

      return(list(
        is_triphasic = TRUE,
        triphasic_type = "exc_inh_exc",
        e1_start = e1$start, e1_end = e1$end,
        i_start  = inh$start, i_end  = inh$end,
        e2_start = e2$start, e2_end = e2$end,
        reason = NA_character_
      ))
    }
  }

  list(
    is_triphasic = FALSE,
    triphasic_type = NA_character_,
    e1_start = NA_real_, e1_end = NA_real_,
    i_start  = NA_real_, i_end  = NA_real_,
    e2_start = NA_real_, e2_end = NA_real_,
    reason = "no_valid_ordering"
  )
}
add_triphasic_features <- function(df_prepped,
                                   per_cell,
                                   window_total = c(0, 20),
                                   thr = 10,
                                   k_smooth = 11,
                                   smooth_method = c("rollmean","rollmedian","loess"),
                                   loess_span = 0.15,
                                   min_persist_s = 0.5,
                                   min_gap_s = 0.05) {

  smooth_method <- match.arg(smooth_method)

  # compute triphasic results per cell
  tri_tbl <- df_prepped %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(function(.x, .y) {
      d <- .x %>%
        dplyr::arrange(.data$time) %>%
        add_smoothed_trace(k_smooth = k_smooth, smooth_method = smooth_method, loess_span = loess_span)

      tri <- detect_triphasic_exc_inh_exc(
        d,
        time_col = "time",
        y_col = "pct_smooth",
        window_total = window_total,
        thr = thr,
        min_persist_s = min_persist_s,
        min_gap_s = min_gap_s,
        baseline = 100
      )

      tibble::tibble(
        is_triphasic = tri$is_triphasic,
        triphasic_type = tri$triphasic_type,
        tri_e1_start = tri$e1_start, tri_e1_end = tri$e1_end,
        tri_i_start  = tri$i_start,  tri_i_end  = tri$i_end,
        tri_e2_start = tri$e2_start, tri_e2_end = tri$e2_end,
        tri_reason = tri$reason,
        tri_thr = thr,
        tri_window0 = window_total[1],
        tri_window1 = window_total[2],
        tri_min_persist_s = min_persist_s,
        tri_min_gap_s = min_gap_s
      )
    }) %>%
    dplyr::ungroup()

  per_cell %>%
    dplyr::left_join(tri_tbl, by = "cell_name")
}
per_cell2 <- add_triphasic_features(
  df_prepped = df_prepped,
  per_cell   = per_cell2,
  window_total = c(0, 60),      # adjust if your rebound is later/earlier
  thr = 10,                     # matches min_excursion_pct scale
  k_smooth = 11,
  smooth_method = "rollmean",
  min_persist_s = 0.5,
  min_gap_s = 0.05
)


per_cell3 <- per_cell2 %>%
  dplyr::mutate(
    response_motif = dplyr::case_when(
      is_triphasic ~ "triphasic_exc_inh_exc",
      response_type == "biphasic" ~ "biphasic",
      response_type %in% c("excitation","inhibition") ~ as.character(response_type),
      TRUE ~ "no_change"
    )
  )

d <- df_prepped %>%
  dplyr::filter(cell_name == CELL_ID) %>%
  dplyr::arrange(time) %>%
  add_smoothed_trace(
    k_smooth = 11,
    smooth_method = "rollmean",
    loess_span = 0.15
  )

# example inside plot_bucket_qc() after you have d with pct_smooth
tri <- detect_triphasic_exc_inh_exc(d, window_total = c(0, 20), thr = min_excursion_pct,
                                    min_persist_s = 0.5, min_gap_s = 0.05)

if (isTRUE(tri$is_triphasic)) {
  phase_df <- tibble::tibble(
    phase = c("E1","I","E2"),
    xmin = c(tri$e1_start, tri$i_start, tri$e2_start),
    xmax = c(tri$e1_end,   tri$i_end,   tri$e2_end)
  )
  p <- p + ggplot2::geom_rect(
    data = phase_df,
    ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE,
    alpha = 0.06
  )
}
