#' #' Plot pulse-train summary features for an RMP protocol
#' #'
#' #' Creates a faceted feature-by-pulse plot from the output of
#' #' \code{analyze_rmp_protocol_set()} (specifically \code{proto_analysis$pulse}).
#' #'
#' #' @param pulse A list like \code{proto_analysis$pulse} containing tables
#' #'   \code{pulse_on}, \code{pulse_tau}, \code{pulse_off}. (May also contain \code{seg_info}.)
#' #' @param title Optional plot title.
#' #' @param subtitle Optional plot subtitle.
#' #' @param point_size Point size for geom_point.
#' #' @param line_width Line width for geom_line.
#' #' @param save_png Logical; if TRUE save to \code{png_path}.
#' #' @param png_path Output path for PNG if \code{save_png=TRUE}.
#' #' @param dpi PNG dpi.
#' #' @param width,height PNG size in inches.
#' #'
#' #' @return A list with:
#' #' \describe{
#' #'   \item{plot}{ggplot object}
#' #'   \item{df_wide}{wide joined feature table}
#' #'   \item{df_long}{long format used for plotting}
#' #' }
#' #' @export
#' plot_rmp_pulse_summary <- function(pulse,
#'                                    title = NULL,
#'                                    subtitle = NULL,
#'                                    point_size = 1.4,
#'                                    line_width = 0.4,
#'                                    save_png = FALSE,
#'                                    png_path = NULL,
#'                                    dpi = 300,
#'                                    width = 7.5,
#'                                    height = 8) {
#'
#'   if (is.null(pulse) || !is.list(pulse)) {
#'     stop("pulse must be a list like proto_analysis$pulse.")
#'   }
#'
#'   pulse_on  <- pulse$pulse_on  %||% NULL
#'   pulse_tau <- pulse$pulse_tau %||% NULL
#'   pulse_off <- pulse$pulse_off %||% NULL
#'
#'   if (is.null(pulse_on) && is.null(pulse_off)) {
#'     stop("pulse appears empty: need pulse_on and/or pulse_off tables.")
#'   }
#'
#'   # ---------- helpers ----------
#'   get_idx <- function(id, prefix) {
#'     # id like "on_003" or "off_012"
#'     out <- suppressWarnings(as.integer(sub(paste0("^", prefix, "_"), "", as.character(id))))
#'     out
#'   }
#'
#'   # ---------- build ON wide ----------
#'   on_wide <- NULL
#'   if (!is.null(pulse_on) && nrow(pulse_on) > 0) {
#'
#'     on_wide <- pulse_on |>
#'       dplyr::mutate(
#'         pulse_state = "pulse_on",
#'         pulse_idx = get_idx(.data$pulse_id, "on")
#'       ) |>
#'       dplyr::select(.data$pulse_state, .data$pulse_idx,
#'                     mean_on = .data$mean_on,
#'                     peak_neg = .data$peak_neg,
#'                     t_peak = .data$t_peak,
#'                     duration = .data$duration,
#'                     delta_mean_minus_peak = .data$delta_mean_minus_peak)
#'
#'     if (!is.null(pulse_tau) && nrow(pulse_tau) > 0) {
#'       tau_wide <- pulse_tau |>
#'         dplyr::mutate(pulse_idx = get_idx(.data$pulse_id, "on")) |>
#'         dplyr::select(.data$pulse_idx, tau = .data$tau)
#'
#'       on_wide <- on_wide |>
#'         dplyr::left_join(tau_wide, by = "pulse_idx")
#'     }
#'   }
#'
#'   # ---------- build OFF wide ----------
#'   off_wide <- NULL
#'   if (!is.null(pulse_off) && nrow(pulse_off) > 0) {
#'
#'     off_wide <- pulse_off |>
#'       dplyr::mutate(
#'         pulse_state = "pulse_off",
#'         pulse_idx = get_idx(.data$pulse_id, "off")
#'       ) |>
#'       dplyr::select(.data$pulse_state, .data$pulse_idx,
#'                     mean_off = .data$mean_off,
#'                     early_peak_off = .data$early_peak_off,
#'                     early_t_peak_off = .data$early_t_peak_off,
#'                     duration = .data$duration)
#'   }
#'
#'   df_wide <- dplyr::bind_rows(on_wide, off_wide)
#'
#'   if (is.null(df_wide) || nrow(df_wide) == 0) {
#'     stop("No rows to plot after processing pulse tables.")
#'   }
#'
#'   # ---------- long format ----------
#'   # Make nice feature labels (matching your figure style)
#'   feature_labels <- c(
#'     peak_neg = "Pulse ON: peak negative (mV)",
#'     mean_on = "Pulse ON: mean (mV)",
#'     delta_mean_minus_peak = "Pulse ON: mean \u2212 peak (mV)",
#'     t_peak = "Pulse ON: time to peak (s)",
#'     tau = "Pulse ON: tau to steady (s)",
#'     mean_off = "Pulse OFF: mean (mV)",
#'     early_peak_off = paste0("Pulse OFF: early peak \u2264", round((pulse_off$early_window_s %||% 0.1)[1] * 1000), "ms (mV)"),
#'     early_t_peak_off = "Pulse OFF: time to early peak (s)"
#'   )
#'
#'   # pivot_longer over numeric feature columns except pulse_idx/state
#'   df_long <- df_wide |>
#'     dplyr::mutate(pulse_idx = as.integer(.data$pulse_idx)) |>
#'     tidyr::pivot_longer(
#'       cols = -c(.data$pulse_state, .data$pulse_idx),
#'       names_to = "feature",
#'       values_to = "value"
#'     ) |>
#'     dplyr::filter(is.finite(.data$pulse_idx)) |>
#'     dplyr::filter(is.finite(.data$value)) |>
#'     dplyr::mutate(
#'       feature_label = dplyr::if_else(
#'         .data$feature %in% names(feature_labels),
#'         unname(feature_labels[.data$feature]),
#'         .data$feature
#'       ),
#'       # order panels: ON block then OFF block
#'       panel_order = dplyr::case_when(
#'         grepl("^Pulse ON", .data$feature_label) ~ 1L,
#'         grepl("^Pulse OFF", .data$feature_label) ~ 2L,
#'         TRUE ~ 3L
#'       )
#'     ) |>
#'     dplyr::arrange(.data$panel_order, .data$feature_label, .data$pulse_idx)
#'
#'   # Keep panel ordering stable
#'   df_long$feature_label <- factor(df_long$feature_label,
#'                                   levels = unique(df_long$feature_label))
#'
#'   # ---------- plot ----------
#'   p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data$pulse_idx, y = .data$value, color = .data$pulse_state)) +
#'     ggplot2::geom_point(size = point_size, alpha = 0.85) +
#'     ggplot2::geom_line(ggplot2::aes(group = .data$pulse_state), linewidth = line_width, alpha = 0.6) +
#'     ggplot2::facet_grid(.data$feature_label ~ ., scales = "free_y", switch = "y") +
#'     ggplot2::theme_classic() +
#'     ggplot2::labs(
#'       x = "Pulse number",
#'       y = NULL,
#'       color = "State",
#'       title = title,
#'       subtitle = subtitle
#'     ) +
#'     ggplot2::theme(
#'       legend.position = "top",
#'       strip.background = ggplot2::element_blank(),
#'       strip.placement = "outside",
#'       strip.text.y.left = ggplot2::element_text(angle = 0)
#'     )
#'
#'   if (isTRUE(save_png)) {
#'     if (is.null(png_path) || !nzchar(png_path)) {
#'       stop("save_png=TRUE requires png_path.")
#'     }
#'     ggplot2::ggsave(filename = png_path, plot = p, dpi = dpi, width = width, height = height, units = "in")
#'   }
#'
#'   list(plot = p, df_wide = df_wide, df_long = df_long)
#' }
