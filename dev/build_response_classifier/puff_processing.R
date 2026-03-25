# puff_processing.R

#source("puff_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

prepare_puff_df <- function(df,
                            baseline_window = c(-5, 0),
                            subclass_keep = c("L23_IT", "L5_ET", "L5_IT"),
                            recode_subclass = c("L3c" = "L23_IT"),
                            bath_keep = "none",
                            puff_pattern = "^5HT") {

  assert_has_cols(
    df,
    c("cell_name", "assigned_subclass", "time", "x_bin", "instRate",
      "bath", "puff")
  )

  df0 <- df %>%
    dplyr::mutate(
      time = as.numeric(.data$time),
      x_bin = as.numeric(.data$x_bin),
      instRate = as.numeric(.data$instRate),
      assigned_subclass = dplyr::if_else(
        .data$assigned_subclass %in% names(recode_subclass),
        unname(recode_subclass[.data$assigned_subclass]),
        .data$assigned_subclass
      )
    ) %>%
    dplyr::filter(
      .data$bath == bath_keep,
      .data$assigned_subclass %in% subclass_keep,
      grepl(puff_pattern, .data$puff)
    )

  base_tbl <- df0 %>%
    dplyr::filter(
      .data$time >= baseline_window[1],
      .data$time < baseline_window[2]
    ) %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::summarise(
      baseline_hz = mean(.data$instRate, na.rm = TRUE),
      n_baseline_pts = sum(is.finite(.data$instRate)),
      .groups = "drop"
    )

  df_prepped <- df0 %>%
    dplyr::left_join(base_tbl, by = "cell_name") %>%
    dplyr::mutate(
      pct_of_baseline = dplyr::if_else(
        is.finite(.data$baseline_hz) & .data$baseline_hz > 0,
        100 * .data$instRate / .data$baseline_hz,
        NA_real_
      )
    )

  df_prepped
}

extract_cell_features_one <- function(df_cell,
                                      rapid_window = c(0, 10),
                                      mid_window = c(10, 40),
                                      late_window = c(40, 60),
                                      post_window = c(0, 60),
                                      baseline = 100,
                                      min_excursion_pct = 10,
                                      min_persist_s = 0.5,
                                      smooth_method = "rollmean",
                                      smooth_k = 3,
                                      loess_span = 0.2) {

  assert_has_cols(df_cell, c("time", "x_bin", "instRate", "pct_of_baseline"))

  df_cell <- df_cell %>%
    dplyr::arrange(.data$time)

  time <- df_cell$time
  y_raw <- df_cell$pct_of_baseline
  y_smooth <- smooth_vec(
    y = y_raw,
    method = smooth_method,
    k = smooth_k,
    loess_span = loess_span,
    x = time
  )

  rapid_min <- get_extreme_in_window(time, y_smooth, rapid_window, "min")
  rapid_max <- get_extreme_in_window(time, y_smooth, rapid_window, "max")
  abs_min   <- get_extreme_in_window(time, y_smooth, post_window, "min")
  abs_max   <- get_extreme_in_window(time, y_smooth, post_window, "max")

  rapid_has_exc <- has_persistent_excursion(
    time, y_smooth,
    window = rapid_window,
    direction = "exc",
    baseline = baseline,
    thr = min_excursion_pct,
    min_persist_s = min_persist_s
  )

  rapid_has_inh <- has_persistent_excursion(
    time, y_smooth,
    window = rapid_window,
    direction = "inh",
    baseline = baseline,
    thr = min_excursion_pct,
    min_persist_s = min_persist_s
  )

  abs_has_exc <- has_persistent_excursion(
    time, y_smooth,
    window = post_window,
    direction = "exc",
    baseline = baseline,
    thr = min_excursion_pct,
    min_persist_s = min_persist_s
  )

  abs_has_inh <- has_persistent_excursion(
    time, y_smooth,
    window = post_window,
    direction = "inh",
    baseline = baseline,
    thr = min_excursion_pct,
    min_persist_s = min_persist_s
  )

  max_is_real <- is.finite(abs_max$value) &&
    abs_max$value >= (baseline + min_excursion_pct) &&
    abs_has_exc

  min_is_real <- is.finite(abs_min$value) &&
    abs_min$value <= (baseline - min_excursion_pct) &&
    abs_has_inh

  exc_mag <- if (max_is_real) abs_max$value - baseline else NA_real_
  inh_mag <- if (min_is_real) baseline - abs_min$value else NA_real_

  direction <- dplyr::case_when(
    is.finite(exc_mag) & is.finite(inh_mag) ~ ifelse(exc_mag >= inh_mag, "excitation", "inhibition"),
    is.finite(exc_mag) ~ "excitation",
    is.finite(inh_mag) ~ "inhibition",
    TRUE ~ "no_change"
  )

  response_peak <- dplyr::case_when(
    direction == "excitation" ~ exc_mag,
    direction == "inhibition" ~ inh_mag,
    TRUE ~ NA_real_
  )

  response_type <- dplyr::case_when(
    rapid_has_exc & rapid_has_inh ~ "biphasic",
    max_is_real & !min_is_real ~ "excitation",
    min_is_real & !max_is_real ~ "inhibition",
    max_is_real & min_is_real ~ "biphasic",
    TRUE ~ "no_change"
  )

  biphasic_order <- dplyr::case_when(
    !(rapid_has_exc & rapid_has_inh) ~ NA_character_,
    !is.finite(rapid_max$time) | !is.finite(rapid_min$time) ~ "ambiguous",
    rapid_max$time < rapid_min$time ~ "exc_then_inh",
    rapid_min$time < rapid_max$time ~ "inh_then_exc",
    TRUE ~ "ambiguous"
  )

  biphasic_first <- dplyr::case_when(
    biphasic_order == "exc_then_inh" ~ "excitation_first",
    biphasic_order == "inh_then_exc" ~ "inhibition_first",
    TRUE ~ NA_character_
  )

  mid_mean <- mean(
    y_smooth[is.finite(time) & time >= mid_window[1] & time <= mid_window[2]],
    na.rm = TRUE
  )
  if (!is.finite(mid_mean)) mid_mean <- NA_real_

  late_mean <- mean(
    y_smooth[is.finite(time) & time >= late_window[1] & time <= late_window[2]],
    na.rm = TRUE
  )
  if (!is.finite(late_mean)) late_mean <- NA_real_

  dlate <- tibble::tibble(time = time, y = y_smooth) %>%
    dplyr::filter(
      is.finite(.data$time),
      is.finite(.data$y),
      .data$time >= late_window[1],
      .data$time <= late_window[2]
    )

  late_slope <- if (nrow(dlate) >= 5) {
    stats::coef(stats::lm(y ~ time, data = dlate))[["time"]]
  } else {
    NA_real_
  }

  mag_for_call <- dplyr::case_when(
    response_type == "excitation" ~ exc_mag,
    response_type == "inhibition" ~ inh_mag,
    TRUE ~ NA_real_
  )

  mag_bin <- classify_magnitude(mag_for_call)

  response_call <- dplyr::case_when(
    response_type == "excitation" & mag_bin == "mild" ~ "mild_excitation",
    response_type == "excitation" & mag_bin == "moderate" ~ "moderate_excitation",
    response_type == "excitation" & mag_bin == "strong" ~ "strong_excitation",
    response_type == "inhibition" & mag_bin == "mild" ~ "mild_inhibition",
    response_type == "inhibition" & mag_bin == "moderate" ~ "moderate_inhibition",
    response_type == "inhibition" & mag_bin == "strong" ~ "strong_inhibition",
    TRUE ~ "no_change"
  )

  tibble::tibble(
    assigned_subclass = first_non_missing(df_cell$assigned_subclass),
    predicted_subclass = first_non_missing(df_cell$predicted_subclass),
    Species = first_non_missing(df_cell$Species),
    blockers = first_non_missing(df_cell$blockers),
    bath = first_non_missing(df_cell$bath),
    puff = first_non_missing(df_cell$puff),
    expCon = first_non_missing(df_cell$expCon),
    assigned_depth = suppressWarnings(as.numeric(first_non_missing(df_cell$assigned_depth))),
    cluster_Corr = first_non_missing(df_cell$cluster_Corr),

    baseline_hz = first_non_missing(df_cell$baseline_hz, default = NA_real_),
    n_timepoints = nrow(df_cell),

    rapid_min_time = rapid_min$time,
    rapid_min_pct = rapid_min$value,
    rapid_max_time = rapid_max$time,
    rapid_max_pct = rapid_max$value,

    abs_min_time = abs_min$time,
    abs_min_pct = abs_min$value,
    abs_max_time = abs_max$time,
    abs_max_pct = abs_max$value,

    rapid_has_exc = rapid_has_exc,
    rapid_has_inh = rapid_has_inh,
    min_is_real = min_is_real,
    max_is_real = max_is_real,

    exc_mag = exc_mag,
    inh_mag = inh_mag,
    direction = direction,
    response_peak = response_peak,

    response_type = response_type,
    biphasic_order = biphasic_order,
    biphasic_first = biphasic_first,
    response_call = response_call,

    mid_mean_pct = mid_mean,
    late_mean_pct = late_mean,
    late_slope_pct_per_s = late_slope,

    smooth_method = smooth_method,
    smooth_k = smooth_k,
    min_excursion_pct = min_excursion_pct,
    min_persist_s = min_persist_s
  )
}

build_per_cell_classifier <- function(df_prepped,
                                      rapid_window = c(0, 10),
                                      mid_window = c(10, 40),
                                      late_window = c(40, 60),
                                      post_window = c(0, 60),
                                      min_excursion_pct = 10,
                                      min_persist_s = 0.5,
                                      smooth_method = "rollmean",
                                      smooth_k = 3,
                                      loess_span = 0.2) {

  assert_has_cols(
    df_prepped,
    c("cell_name", "time", "x_bin", "instRate", "pct_of_baseline", "baseline_hz")
  )

  out <- df_prepped %>%
    dplyr::group_by(.data$cell_name) %>%
    dplyr::group_modify(function(.x, .y) {
      extract_cell_features_one(
        df_cell = .x,
        rapid_window = rapid_window,
        mid_window = mid_window,
        late_window = late_window,
        post_window = post_window,
        min_excursion_pct = min_excursion_pct,
        min_persist_s = min_persist_s,
        smooth_method = smooth_method,
        smooth_k = smooth_k,
        loess_span = loess_span
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      response_call = factor(
        .data$response_call,
        levels = c(
          "strong_inhibition",
          "moderate_inhibition",
          "mild_inhibition",
          "no_change",
          "mild_excitation",
          "moderate_excitation",
          "strong_excitation"
        )
      ),
      response_type = factor(
        .data$response_type,
        levels = c("inhibition", "excitation", "biphasic", "no_change")
      )
    )

  out
}
