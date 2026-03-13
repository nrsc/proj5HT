## =========================================================
## TRACE-SAFE CLONE of your reference build (minimal diffs)
## - baseline computed per trace_id
## - feature extraction grouped per trace_id
## - NO extra baseline gating beyond baseline_hz > 0
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(stringr)
})

# ---- make_trace_id: same idea you’ve already accepted ----
make_trace_id <- function(df, key_cols = c("cell_name","puff","bath","blockers","protocol")) {
  use <- key_cols[key_cols %in% names(df)]
  if (!("cell_name" %in% use)) use <- c("cell_name", use)

  df %>%
    tidyr::unite("trace_id", dplyr::all_of(use), sep = "|", remove = FALSE, na.rm = TRUE) %>%
    dplyr::pull(trace_id)
}

# =========================================================
# 1) Prep raw df: baseline and pct_of_baseline (PER TRACE)
#    Faithful rule: baseline_hz > 0 only
# =========================================================
prep_puff_df_trace <- function(df,
                               baseline_window = c(-5, 0)) {

  stopifnot(all(c("cell_name","time","instRate") %in% names(df)))

  df0 <- df %>%
    dplyr::mutate(
      time = as.numeric(time),
      instRate = as.numeric(instRate)
    )

  # ensure optional cols exist for downstream plotting
  if (!("blockers" %in% names(df0))) df0 <- df0 %>% dplyr::mutate(blockers = NA_character_)
  if (!("assigned_subclass" %in% names(df0))) df0 <- df0 %>% dplyr::mutate(assigned_subclass = NA_character_)

  # add trace_id
  df0 <- df0 %>% dplyr::mutate(trace_id = make_trace_id(df0))

  base_tbl <- df0 %>%
    dplyr::filter(time >= baseline_window[1], time < baseline_window[2]) %>%
    dplyr::group_by(trace_id) %>%
    dplyr::summarise(
      baseline_hz = mean(instRate, na.rm = TRUE),
      n_baseline_pts = sum(is.finite(instRate)),
      .groups = "drop"
    )

  df1 <- df0 %>%
    dplyr::left_join(base_tbl, by = "trace_id") %>%
    dplyr::mutate(
      pct_of_baseline = dplyr::if_else(
        is.finite(baseline_hz) & baseline_hz > 0,
        100 * instRate / baseline_hz,
        NA_real_
      )
    )

  df1
}

# =========================================================
# 2) Smoothing helper (UNCHANGED except explicit dplyr::)
# =========================================================
add_smoothed_trace <- function(df_cell,
                               k_smooth = 11,
                               smooth_method = c("rollmedian","rollmean","loess"),
                               loess_span = 0.15) {
  smooth_method <- match.arg(smooth_method)

  df_cell <- df_cell %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(time = as.numeric(time))

  y <- df_cell$pct_of_baseline

  y_smooth <- switch(
    smooth_method,
    rollmedian = zoo::rollmedian(y, k = k_smooth, fill = NA, align = "center"),
    rollmean   = zoo::rollmean(y,   k = k_smooth, fill = NA, align = "center"),
    loess = {
      ok <- is.finite(df_cell$time) & is.finite(y)
      if (sum(ok) < 10) {
        rep(NA_real_, length(y))
      } else {
        fit <- stats::loess(y[ok] ~ df_cell$time[ok], span = loess_span, degree = 1)
        out <- rep(NA_real_, length(y))
        out[ok] <- stats::predict(fit, newdata = df_cell$time[ok])
        out
      }
    }
  )

  df_cell %>% dplyr::mutate(pct_smooth = as.numeric(y_smooth))
}

# =========================================================
# 3) has_persistent_excursion (UNCHANGED except explicit dplyr::)
# =========================================================
has_persistent_excursion <- function(d, window, direction = c("exc","inh"),
                                     thr = 5,
                                     col = "pct_smooth",
                                     min_persist_s = 0.5) {
  direction <- match.arg(direction)

  d0 <- d %>%
    dplyr::filter(time >= window[1], time <= window[2]) %>%
    dplyr::arrange(time) %>%
    dplyr::filter(is.finite(.data[[col]]), is.finite(time))

  if (nrow(d0) < 2) return(FALSE)

  flag <- if (direction == "exc") d0[[col]] >= (100 + thr) else d0[[col]] <= (100 - thr)
  if (!any(flag)) return(FALSE)

  r <- rle(flag)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1

  true_runs <- which(r$values)
  if (length(true_runs) == 0) return(FALSE)

  for (k in true_runs) {
    i0 <- starts[k]
    i1 <- ends[k]
    dur <- d0$time[i1] - d0$time[i0]
    if (is.finite(dur) && dur >= min_persist_s) return(TRUE)
  }

  FALSE
}

# =========================================================
# 4) extract_bucket_features_one (KEEP YOUR EXACT VERSION)
#    IMPORTANT: it uses df_cell$cell_name[1] as the ID
# =========================================================
# --> use your existing extract_bucket_features_one() unchanged

# =========================================================
# 4b) TRACE-SAFE bulk extractor: group by trace_id
#     and temporarily set cell_name = trace_id for extractor
# =========================================================
extract_bucket_features_trace <- function(df_prepped,
                                          rapid_window = c(0, 10),
                                          mid_window   = c(10, 40),
                                          late_slope_window = c(36, 56),
                                          late_mean_window  = c(46, 56),
                                          post_window = c(0, 56),
                                          min_excursion_pct = 10,
                                          k_smooth = 11,
                                          smooth_method = c("rollmedian","rollmean","loess"),
                                          loess_span = 0.15,
                                          min_persist_s = 0.5) {

  smooth_method <- match.arg(smooth_method)

  df_prepped %>%
    dplyr::group_by(trace_id) %>%
    dplyr::group_modify(~{
      d <- .x
      tid <- .y$trace_id[[1]]

      # extractor wants a cell_name identifier; give it trace_id
      d_for <- d %>% dplyr::mutate(cell_name = tid)

      feat <- extract_bucket_features_one(
        d_for,
        rapid_window = rapid_window,
        mid_window = mid_window,
        late_slope_window = late_slope_window,
        late_mean_window = late_mean_window,
        post_window = post_window,
        min_excursion_pct = min_excursion_pct,
        k_smooth = k_smooth,
        smooth_method = smooth_method,
        loess_span = loess_span,
        min_persist_s = min_persist_s
      )

      # return WITHOUT trace_id column (group_modify will add it)
      # keep real cell_name as metadata
      feat %>%
        dplyr::select(-cell_name) %>%
        dplyr::mutate(
          cell_name = d$cell_name[1],
          assigned_subclass = d$assigned_subclass[1],
          blockers = d$blockers[1]
        )
    }) %>%
    dplyr::ungroup()
}

# =========================================================
# 5) add_biphasic_features (USE YOUR EXACT VERSION)
# =========================================================
# --> use your existing add_biphasic_features() unchanged
add_biphasic_features <- function(per_cell) {

  need <- c("rapid_pattern","rapid_min_time","rapid_max_time","rapid_min_pct","rapid_max_pct",
            "abs_min_pct","abs_max_pct")
  miss <- setdiff(need, names(per_cell))
  if (length(miss) > 0) stop("per_cell missing: ", paste(miss, collapse = ", "))

  per_cell %>%
    dplyr::mutate(
      is_biphasic_rapid = rapid_pattern %in% c("rapid_biphasic_exc_to_inh","rapid_biphasic_inh_to_exc"),

      biphasic_order = dplyr::case_when(
        rapid_pattern == "rapid_biphasic_exc_to_inh" ~ "exc_then_inh",
        rapid_pattern == "rapid_biphasic_inh_to_exc" ~ "inh_then_exc",
        TRUE ~ NA_character_
      ),

      biphasic_first = dplyr::case_when(
        biphasic_order == "exc_then_inh" ~ "excitation_first",
        biphasic_order == "inh_then_exc" ~ "inhibition_first",
        TRUE ~ NA_character_
      ),

      # response_type: biphasic takes precedence; then direction
      response_type = dplyr::case_when(
        is_biphasic_rapid ~ "biphasic",
        direction %in% c("excitation","inhibition") ~ direction,
        TRUE ~ "no_change"
      )
    )
}

## ---- user-side filtering (same as you) ----
df_use <- df %>%
  dplyr::filter(
    bath == "none",
    !cell_name %in% c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03"),
    grepl("^5HT", puff)
  ) %>%
  dplyr::mutate(
    assigned_subclass = dplyr::case_when(
      assigned_subclass %in% c("L3c") ~ "L23_IT",
      TRUE ~ assigned_subclass
    )
  ) %>%
  dplyr::filter(assigned_subclass %in% c("L23_IT","L5_ET","L5_IT"))

## ---- prep with baseline [-5,0) ----
df_prepped <- prep_puff_df_trace(df_use, baseline_window = c(-5, 0))

## ---- per-trace features (SMOOTH+PERSISTENCE) ----
per_cell <- extract_bucket_features_trace(
  df_prepped,
  rapid_window = c(0, 10),
  mid_window = c(10, 40),
  late_slope_window = c(36, 56),
  late_mean_window  = c(46, 56),
  post_window = c(0, 56),
  min_excursion_pct = 10,
  k_smooth = 3,
  smooth_method = "rollmean",
  min_persist_s = 0.5
)

## ---- add biphasic labels ----
per_cell2 <- per_cell %>%
  add_biphasic_features() %>%
  dplyr::mutate(
    biphasic_order = factor(biphasic_order,
                            levels = c("exc_then_inh","inh_then_exc")),
    response_type = factor(response_type,
                           levels = c("excitation","inhibition","biphasic","no_change"))
  )

## ---- blockers_bin directly (already present, no join ambiguity) ----
per_cell2 <- per_cell2 %>%
  dplyr::mutate(
    blockers_bin = dplyr::if_else(is.na(blockers) | blockers == "",
                                  "no_blockers", "has_blockers"),
    blockers_bin = factor(blockers_bin, levels = c("no_blockers","has_blockers"))
  )

## sanity:
per_cell2 %>% dplyr::count(response_type, sort = TRUE)
per_cell2 %>% dplyr::filter(response_type=="biphasic") %>% dplyr::count(biphasic_order, sort = TRUE)

