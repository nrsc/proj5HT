#' Plot baseline vs wash ASP protocols side-by-side (trace + heat strip)
#'
#' Builds two single-protocol ASP plots (baseline and wash) using
#' `plot_asp_single_protocol()` and combines them into a two-column layout
#' with consistent styling.
#'
#' The function expects protocol keys to be stored in `x$dfs$dwin$baseline_key`
#' and `x$dfs$dwin$wash_key`. It also expects the ASP object at `x$dfs$asp`.
#'
#' @param x A nested list (e.g., your nnest/srt-like object) that contains:
#'   \itemize{
#'     \item `x$dfs$asp`: ASP object compatible with `plot_asp_single_protocol()`
#'     \item `x$dfs$dwin$baseline_key`: protocol key for baseline
#'     \item `x$dfs$dwin$wash_key`: protocol key for wash
#'   }
#' @param baseline_key Character scalar. Baseline protocol key. Default:
#'   `x$dfs$dwin$baseline_key`.
#' @param wash_key Character scalar. Wash protocol key. Default:
#'   `x$dfs$dwin$wash_key`.
#' @param baseline_heat_dt Numeric scalar. Bin width (s) for baseline heat strip.
#'   Default `0.5`.
#' @param wash_heat_dt Numeric scalar. Bin width (s) for wash heat strip.
#'   Default `1`.
#' @param ttl_shape Numeric/integer. Shape code passed through to
#'   `plot_asp_single_protocol()` for TTL markers. Default `6`.
#' @param y_mode Character. Y scaling mode passed to `plot_asp_single_protocol()`.
#'   Default `"log2_fc"`.
#' @param y_lim Numeric length-2 vector. Shared y-limits for both panels.
#'   Default `c(-2, 1)`.
#' @param save_png Logical. Passed to `plot_asp_single_protocol()` for both plots.
#'   Default `FALSE`.
#' @param show_plot Logical. Passed to `plot_asp_single_protocol()` for both plots.
#'   Default `FALSE`.
#'
#' @return A `ggplot`/`patchwork` object containing the two plots side-by-side.
#'
#' @examples
#' \dontrun{
#' p <- plot_sidebyside_heat_wash(x = srt)
#' p
#'
#' # Override defaults
#' p2 <- plot_sidebyside_heat_wash(
#'   x = srt,
#'   heat_dt = 1,
#'   heat_dt = 1,
#'   y_lim = c(-3, 2)
#' )
#' }
#'
#' @export
plot_sidebyside_heat_wash <- function(
    x,
    baseline_key = x$dfs$dwin$baseline_key,
    wash_key     = x$dfs$dwin$wash_key,
    heat_dt = 1,
    ttl_shape = 6,
    y_mode = "log2_fc",
    y_lim  = c(-2, 1),
    save_png  = FALSE,
    show_plot = FALSE
){

  resBl <- plot_asp_single_protocol(
    x$dfs$asp,
    protocol_key = baseline_key,
    key_match = TRUE,
    heat_dt = heat_dt,
    ttl_shape = ttl_shape,
    y_mode = y_mode,
    y_lim = y_lim,
    save_png = save_png,
    show_plot = show_plot
  )

  if(is.null(resBl)){
    return(NULL)
  }

  resWash <- plot_asp_single_protocol(
    x$dfs$asp,
    protocol_key = wash_key,
    key_match = TRUE,
    heat_dt = heat_dt,
    ttl_shape = ttl_shape,
    y_mode = y_mode,
    y_lim = y_lim,
    save_png = save_png,
    show_plot = show_plot
  )

  p <- (resBl$plot | resWash$plot) +
    patchwork::plot_layout(widths = c(1, 1)) &
    ggplot2::theme(
      plot.subtitle = ggplot2::element_blank(),   # ensure no subtitles anywhere
      plot.title = ggplot2::element_text(size = 11, face = "bold"),
      plot.title.position = "plot"
    )

  return(p)
}
