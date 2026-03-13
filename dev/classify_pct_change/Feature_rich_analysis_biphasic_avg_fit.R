suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# Helpers
# -----------------------------
smooth_ma <- function(x, k = 7) {
  as.numeric(stats::filter(x, rep(1 / k, k), sides = 2))
}

linfit_stats <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 3) return(list(slope = NA_real_, r2 = NA_real_, p = NA_real_))
  fit <- stats::lm(y ~ x)
  s <- summary(fit)
  list(
    slope = unname(coef(fit)[2]),
    r2    = s$r.squared,
    p     = s$coefficients[2, 4]
  )
}

# -----------------------------
# Compute baseline and % of baseline
# -----------------------------
compute_pc <- function(df_cell,
                       baseline_window = c(-5, 0),
                       min_baseline_hz = 0.2,
                       k_smooth = 7) {

  df_cell <- df_cell %>%
    mutate(time = as.numeric(time)) %>%
    arrange(time)

  base <- df_cell %>%
    dplyr::filter(.data$time >= baseline_window[1], .data$time <= baseline_window[2]) %>%
    summarise(b = mean(.data$instRate, na.rm = TRUE)) %>%
    pull(b)

  if (is.na(base) || base < min_baseline_hz) {
    return(list(valid = FALSE, reason = "baseline too low/NA", baseline_hz = base, df = df_cell))
  }

  d <- df_cell %>%
    mutate(
      baseline_hz = base,
      pc_raw    = 100 * .data$instRate / base,
      pc_smooth = smooth_ma(.data$pc_raw, k = k_smooth)
    )

  list(valid = TRUE, reason = NA_character_, baseline_hz = base, df = d)
}

# -----------------------------
# Bucketed response metrics (robust window coercion)
# -----------------------------
bucket_metrics <- function(d,
                           early_window = c(0, 10),
                           mid_window   = c(10, 40),
                           late_window  = c(40, 50),
                           ignore_after_puff = 0.25,
                           baseline_level = 100,
                           min_excursion_pct = 5,
                           summary_fun = c("median", "mean")) {

  # Defensive coercion: fixes your "comparison (<=)..." error
  ew <- as.numeric(unlist(early_window))
  mw <- as.numeric(unlist(mid_window))
  lw <- as.numeric(unlist(late_window))
  if (length(ew) != 2 || length(mw) != 2 || length(lw) != 2) {
    stop("early_window/mid_window/late_window must each be length-2 numeric vectors.")
  }

  summary_fun <- match.arg(summary_fun)
  summ <- function(x) if (summary_fun == "median") stats::median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)

  valcol <- if ("pc_smooth" %in% names(d)) "pc_smooth" else "pc_raw"

  d_early <- d %>% dplyr::filter(.data$time >= (ew[1] + ignore_after_puff), .data$time <= ew[2])
  d_mid   <- d %>% dplyr::filter(.data$time >= mw[1], .data$time <= mw[2])
  d_late  <- d %>% dplyr::filter(.data$time >= lw[1], .data$time <= lw[2])

  early_min <- d_early %>% slice(which.min(.data[[valcol]]))
  early_max <- d_early %>% slice(which.max(.data[[valcol]]))

  early_min_val <- early_min[[valcol]][1]
  early_max_val <- early_max[[valcol]][1]

  exc_present <- is.finite(early_max_val) && early_max_val >= (baseline_level + min_excursion_pct)
  inh_present <- is.finite(early_min_val) && early_min_val <= (baseline_level - min_excursion_pct)

  exc_mag <- if (exc_present) (early_max_val - baseline_level) else NA_real_
  inh_mag <- if (inh_present) (baseline_level - early_min_val) else NA_real_

  immediate_direction <- dplyr::case_when(
    !is.na(exc_mag) & (is.na(inh_mag) | exc_mag > inh_mag) ~ "excitation",
    !is.na(inh_mag) & (is.na(exc_mag) | inh_mag >= exc_mag) ~ "inhibition",
    TRUE ~ "no_change"
  )

  immediate_peak <- dplyr::case_when(
    immediate_direction == "excitation" ~ early_max_val,
    immediate_direction == "inhibition" ~ early_min_val,
    TRUE ~ baseline_level
  )

  immediate_peak_time <- dplyr::case_when(
    immediate_direction == "excitation" ~ early_max$time[1],
    immediate_direction == "inhibition" ~ early_min$time[1],
    TRUE ~ NA_real_
  )

  mid_level <- if (nrow(d_mid) >= 3) summ(d_mid[[valcol]]) else NA_real_
  mid_fit   <- if (nrow(d_mid) >= 5) linfit_stats(d_mid$time, d_mid[[valcol]]) else list(slope = NA_real_, r2 = NA_real_, p = NA_real_)

  mid_direction <- dplyr::case_when(
    is.finite(mid_level) & mid_level >= (baseline_level + min_excursion_pct) ~ "excitation",
    is.finite(mid_level) & mid_level <= (baseline_level - min_excursion_pct) ~ "inhibition",
    is.finite(mid_level) ~ "no_change",
    TRUE ~ NA_character_
  )

  late_mean <- if (nrow(d_late) >= 3) mean(d_late[[valcol]], na.rm = TRUE) else NA_real_
  late_fit  <- if (nrow(d_late) >= 5) linfit_stats(d_late$time, d_late[[valcol]]) else list(slope = NA_real_, r2 = NA_real_, p = NA_real_)

  list(
    early = list(
      min_time = early_min$time[1], min_val = early_min_val, min_flag = if (!inh_present) "min_above_inhibition_threshold" else NA_character_,
      max_time = early_max$time[1], max_val = early_max_val, max_flag = if (!exc_present) "max_below_excursion_threshold" else NA_character_,
      immediate_direction = immediate_direction,
      immediate_peak = immediate_peak,
      immediate_peak_time = immediate_peak_time
    ),
    mid = list(level = mid_level, direction = mid_direction, slope = mid_fit$slope, r2 = mid_fit$r2, p = mid_fit$p),
    late = list(mean = late_mean, slope = late_fit$slope, r2 = late_fit$r2, p = late_fit$p),
    windows = list(early = ew, mid = mw, late = lw)
  )
}

# -----------------------------
# Plot: bucket shading + peak markers + late fit line
# -----------------------------
plot_bucket_qc <- function(df,
                           cell_id,
                           baseline_window = c(-5, 0),
                           early_window = c(0, 10),
                           mid_window   = c(10, 40),
                           late_window  = c(40, 50),
                           ignore_after_puff = 0.25,
                           k_smooth = 7,
                           min_excursion_pct = 5,
                           summary_fun = "median",
                           drop_protocol_baseline = FALSE) {

  df_cell <- df %>% dplyr::filter(.data$cell_name == cell_id)

  if (drop_protocol_baseline && "protocol" %in% names(df_cell)) {
    df_cell <- df_cell %>% dplyr::filter(tolower(.data$protocol) != "baseline")
  }

  out <- compute_pc(df_cell, baseline_window = baseline_window, k_smooth = k_smooth)
  if (!isTRUE(out$valid)) stop(paste("Invalid:", out$reason))

  d <- out$df
  m <- bucket_metrics(
    d,
    early_window = early_window,
    mid_window   = mid_window,
    late_window  = late_window,
    ignore_after_puff = ignore_after_puff,
    min_excursion_pct = min_excursion_pct,
    summary_fun = summary_fun
  )

  # COLORS (as requested): red=excitation, blue=inhibition
  col_exc <- "#B2182B"  # red
  col_inh <- "#2166AC"  # blue

  # Plot peak dots only if threshold passed
  pt_inh <- if (is.na(m$early$min_flag)) data.frame(time = m$early$min_time, y = m$early$min_val) else NULL
  pt_exc <- if (is.na(m$early$max_flag)) data.frame(time = m$early$max_time, y = m$early$max_val) else NULL

  # Late fit line
  lw <- m$windows$late
  d_late <- d %>% dplyr::filter(.data$time >= lw[1], .data$time <= lw[2])
  fit_line <- NULL
  if (nrow(d_late) >= 5) {
    fit <- lm(pc_smooth ~ time, data = d_late)
    fit_line <- data.frame(time = seq(min(d_late$time), max(d_late$time), length.out = 60))
    fit_line$pc_smooth <- predict(fit, newdata = fit_line)
  }

  sub <- paste0(
    "baseline=", round(out$baseline_hz, 2), " Hz | ",
    "early: ", m$early$immediate_direction, " peak=", round(m$early$immediate_peak, 1),
    "% @ ", round(m$early$immediate_peak_time, 2), "s",
    " | mid: ", round(m$mid$level, 1), "% (", m$mid$direction, ")",
    " | late: mean=", round(m$late$mean, 1), "% slope=", signif(m$late$slope, 3), " r2=", signif(m$late$r2, 3)
  )

  ggplot(d, aes(.data$time)) +
    annotate("rect", xmin = baseline_window[1], xmax = baseline_window[2], ymin = -Inf, ymax = Inf,
             alpha = 0.08, fill = "grey60") +
    annotate("rect", xmin = early_window[1], xmax = early_window[2], ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = col_exc) +
    annotate("rect", xmin = mid_window[1], xmax = mid_window[2], ymin = -Inf, ymax = Inf,
             alpha = 0.04, fill = "grey50") +
    annotate("rect", xmin = late_window[1], xmax = late_window[2], ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = col_inh) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40") +
    geom_line(aes(y = .data$pc_raw), color = "grey80", linewidth = 0.4) +
    geom_line(aes(y = .data$pc_smooth), color = "black", linewidth = 0.9) +
    {if (!is.null(pt_inh)) geom_point(data = pt_inh, aes(x = time, y = y), color = col_inh, size = 3)} +
    {if (!is.null(pt_exc)) geom_point(data = pt_exc, aes(x = time, y = y), color = col_exc, size = 3)} +
    {if (!is.null(fit_line)) geom_line(data = fit_line, aes(x = time, y = pc_smooth), linewidth = 1.2)} +
    labs(
      title = paste0(cell_id, " | ", unique(d$assigned_subclass)[1]),
      subtitle = sub,
      x = "Time (s)",
      y = "% of baseline (100 = baseline)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 8)
    )
}

#Example
p <- plot_bucket_qc(df, "QF25.26.024.19.07.04", k_smooth = 11, min_excursion_pct = 5, ignore_after_puff = 0.5)
print(p)
p <- plot_bucket_qc(df, "QF25.26.022.19.09.01", k_smooth = 40, min_excursion_pct = 5)

save_bucket_qc_batch <- function(df, cell_ids, out_dir = "qc_bucketed",
                                 k_smooth = 7, min_excursion_pct = 5) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  for (cid in cell_ids) {
    p <- try(plot_bucket_qc(df, cid, k_smooth = k_smooth, min_excursion_pct = min_excursion_pct), silent = TRUE)
    if (!inherits(p, "try-error")) {
      ggsave(file.path(out_dir, paste0(cid, ".png")), p, width = 7, height = 5, dpi = 150)
    }
  }
}
