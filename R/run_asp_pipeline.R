#' Run the per-cell ASP + wash-in pipeline
#'
#' One-call wrapper around [ensure_asp()] (build / load the Aligned Spike
#' Protocol payload) and [DrugWashIn_from_asp()] (derive ket/WAY wash-in
#' kinetics from the ASP). Optionally produces concatenated and per-protocol
#' ASP plots.
#'
#' Driven by `dev/startup.sh` (alias `update_5ht`) to incrementally analyze
#' any cell in the `Master_UI_data.ods` Master_Macaque sheet that does not yet
#' have a `<cell>-srt.rds` in the rookery.
#'
#' @param cell Character. proj5HT cell identifier (e.g. `"QN26.26.025.20.10.01"`).
#' @param mMD Optional data frame. Pre-loaded Master_UI sheet rows. If `NULL`,
#'   `mMD_path` + `mMD_sheet` are read on demand.
#' @param mMD_path Path to the Master_UI ODS file.
#' @param mMD_sheet Sheet name within `mMD_path`. Defaults to `"Master_Macaque"`.
#' @param ttl_default Default trial-time-length (seconds) when not specified
#'   in the Master_UI row.
#' @param outlier_sd SD threshold for ASP outlier rejection.
#' @param bl_window Numeric length-2. Baseline window (s) relative to puff
#'   onset.
#' @param save_srt Logical. If `TRUE`, persists the child nnest to
#'   `<rookery>/<cell>/<cell>-srt.rds`.
#' @param make_concat_plot Logical. Render the concatenated ASP plot.
#' @param make_single_plots Logical. Render one plot per ASP protocol.
#' @param make_washin Logical. Run [DrugWashIn_from_asp()] and persist the
#'   wash-in CSV.
#' @param verbose Logical. Echo progress from downstream calls.
#'
#' @return An invisible list with elements:
#'   \describe{
#'     \item{`asp`}{The ASP payload list returned by [ensure_asp()].}
#'     \item{`concat`}{ggplot object (only if `make_concat_plot = TRUE`).}
#'     \item{`single`}{Named list of ggplots, one per protocol key
#'       (only if `make_single_plots = TRUE`).}
#'     \item{`washin`}{Wash-in result list (only if `make_washin = TRUE`).}
#'   }
#'
#' @examples
#' \dontrun{
#' run_asp_pipeline("QN26.26.025.20.10.01")
#' }
#'
#' @export
run_asp_pipeline <- function(cell,
                             mMD = NULL,
                             mMD_path = "~/proj5HT/den/serotonin/Data/Master_UI_data.ods",
                             mMD_sheet = "Master_Macaque",
                             ttl_default = 5,
                             outlier_sd = 5,
                             bl_window = c(-5, 0),
                             save_srt = TRUE,
                             make_concat_plot = TRUE,
                             make_single_plots = TRUE,
                             make_washin = TRUE,
                             verbose = TRUE) {

  got <- ensure_asp(
    cell,
    mMD = mMD,
    mMD_path = mMD_path,
    mMD_sheet = mMD_sheet,
    ttl_default = ttl_default,
    outlier_sd = outlier_sd,
    bl_window = bl_window,
    save_srt = save_srt
  )
  srt <- got$srt
  asp <- got$asp

  res <- list(asp = asp)

  if (isTRUE(make_concat_plot)) {
    res$concat <- plot_asp_concatenated(asp, verbose = verbose)
  }

  if (isTRUE(make_single_plots)) {
    keys <- names(asp$by_protocol)
    res$single <- lapply(keys, function(k) {
      plot_asp_single_protocol(asp, protocol_key = k, verbose = verbose)
    })
    names(res$single) <- keys
  }

  if (isTRUE(make_washin)) {
    res$washin <- DrugWashIn_from_asp(srt, mMD = mMD, plot_it = TRUE, write_csv = TRUE, verbose = verbose)
  }

  invisible(res)
}
