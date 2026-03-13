run_asp_ana_pipeline <- function(srt, verbose = TRUE) {
  stop_if_not_srt(srt)

  if (is.null(srt$ana) || !is.list(srt$ana)) srt$ana <- list()
  if (is.null(srt$ana$asp) || !is.list(srt$ana$asp)) srt$ana$asp <- list()

  # Example: one “canonical” asp summary table
  if (!is.null(srt$dfs$asp$smry)) {
    srt$ana$asp$protocol_summary <- srt$dfs$asp$smry
  }

  # Example: a combined “single-protocol” trace table if you want it
  # (This assumes spikePuff_from_asp returns spike_puff_output; adjust to your return schema.)
  sp <- srt$dfs$single$spikePuff
  if (!is.null(sp) && is.list(sp) && !is.null(sp$spike_puff_output)) {
    srt$ana$asp$spikePuff_trace <- sp$spike_puff_output
  }

  # Example: wash-in trace
  wi <- srt$dfs$drug_wash$washin
  if (!is.null(wi) && is.list(wi) && !is.null(wi$spike_puff_output)) {
    srt$ana$asp$washin_trace <- wi$spike_puff_output
  }

  if (verbose) message("ANA pipeline: built ana$asp tables.")

  srt$ana
}
