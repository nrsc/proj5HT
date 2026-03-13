reticulate::source_python('~/proj5HT/py/plotting_pdf.py')

plot_sweeps_py <- function(nwb_path,
                           sweeps,                # c(90) OR c("data_00089_AD0","data_00090_AD0")
                           fmt = c("pdf","png"),
                           out_dir = file.path("figs","traces"),
                           dpi = 300,
                           width = 7, height = 2.2,
                           lw = 0.8,
                           tlim = NULL,

                           # time handling
                           re_zero = TRUE,
                           gap_s = 0,

                           # TTL
                           ttl_time = NULL,        # seconds after stim sweep start
                           ttl_on_sweep_index = 1, # for baseline+stim, stim is index 2 in R but 1 in python 0-based => 1

                           # publication
                           pub_void_axes = FALSE,
                           scalebar_x = NULL,
                           scalebar_y = NULL,
                           scalebar_label_x = NULL,
                           scalebar_label_y = NULL,
                           scalebar_loc = "lower right",

                           # units
                           y_unit_target = NULL,   # e.g. "mV"
                           title = NULL) {

  fmt <- match.arg(fmt)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(nwb_path)) stop("nwb_path does not exist: ", nwb_path)

  sweeps_py <- if (is.character(sweeps)) as.list(sweeps) else as.list(as.integer(sweeps))
  tlim_py   <- if (is.null(tlim)) reticulate::py_none() else as.numeric(tlim)
  title_py  <- if (is.null(title)) reticulate::py_none() else as.character(title)

  ttl_time_py <- if (is.null(ttl_time)) reticulate::py_none() else as.numeric(ttl_time)
  y_unit_py   <- if (is.null(y_unit_target)) reticulate::py_none() else as.character(y_unit_target)

  out_path <- reticulate::py$plot_nwb_sweeps(
    nwb_path = normalizePath(nwb_path),
    sweeps = sweeps_py,
    out_dir = normalizePath(out_dir, mustWork = FALSE),
    fmt = fmt,
    dpi = as.integer(dpi),
    width = as.numeric(width),
    height = as.numeric(height),
    lw = as.numeric(lw),
    tlim = tlim_py,
    title = title_py,

    re_zero = isTRUE(re_zero),
    gap_s = as.numeric(gap_s),

    ttl_time = ttl_time_py,
    ttl_on_sweep_index = as.integer(ttl_on_sweep_index),

    pub_void_axes = isTRUE(pub_void_axes),
    scalebar_x = if (is.null(scalebar_x)) reticulate::py_none() else as.numeric(scalebar_x),
    scalebar_y = if (is.null(scalebar_y)) reticulate::py_none() else as.numeric(scalebar_y),
    scalebar_label_x = if (is.null(scalebar_label_x)) reticulate::py_none() else as.character(scalebar_label_x),
    scalebar_label_y = if (is.null(scalebar_label_y)) reticulate::py_none() else as.character(scalebar_label_y),
    scalebar_loc = as.character(scalebar_loc),

    y_unit_target = y_unit_py
  )

  out_path
}



sapply(per_cell2$cell_name, function(x){

srt = load5HT(x, rehydrate = FALSE)
#plotBasicPuff(srt, xlim = c(-40,50))

if("baseline" %in% names(srt$dfs$spikeTTL$sweep_names)){
  sweeps = c(
    srt$dfs$spikeTTL$sweep_names[["baseline"]],
    srt$dfs$spikeTTL$sweep_names[["spikeTTL"]]
  )
  ttli = 1
}else{
  sweeps = c(
    srt$dfs$spikeTTL$sweep_names[["spikeTTL"]]
  )
  ttli = 0
}

out <- plot_sweeps_py(
  nwb_path = srt$files$nwb,
  sweeps = sweeps,
  fmt = "pdf",
  title = paste(srt$cell, unique(srt$dfs$spikeTTL$spike_puff_output$assigned_subclass), sep = " -- "),
  re_zero = TRUE,
  gap_s = 0,

  ttl_time = srt$dfs$spikeTTL$ttlTime,
  ttl_on_sweep_index = ttli,

  pub_void_axes = TRUE,
  y_unit_target = "mV",
  scalebar_x = 10,
  scalebar_label_x = "10 s",
  scalebar_y = 20,
  scalebar_label_y = "20 mV"
)

message("Saved: ", out)

})

#rstudioapi::viewer(out)

