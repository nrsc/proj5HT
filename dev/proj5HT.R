### Project looper
library(dplyr)
library(tidyr)
library(nphys)
library(projHCT)
library(proj5HT)

#### Setup ####
MD = projHCT::sheets$MD

# Helper functions
source("~/proj5HT/R/spikePuff_from_asp.R", echo = TRUE)
#file.edit("~/proj5HT/R/spikePuff_from_asp.R", echo = TRUE)

source("~/proj5HT/R/DrugWashIn_from_asp.R", echo = TRUE)
#file.edit("~/proj5HT/R/DrugWashIn_from_asp.R", echo = TRUE)

source("~/proj5HT/R/build_asp.R", echo = TRUE)
#file.edit("~/proj5HT/R/build_asp.R", echo = TRUE)

source("~/proj5HT/R/proj5HT_loop_helper_functions.R", echo = TRUE)
#file.edit("~/proj5HT/R/proj5HT_loop_helper_functions.R", echo = TRUE)

source("~/proj5HT/R/rmp_helper_functions.R")
#file.edit("~/proj5HT/R/rmp_helper_functions.R", echo = TRUE)
#reticulate::source_python('~/proj5HT/py/get_epochs_from_stim.py')


#sheets = projHCT::sheets$files$nnest_files
sh5HT = sheets5HT()
#sh5HT = proj5HT::sh5HT
#mstr = sh5HT$mstrMac
#cells = unique(mstr$cell_name)


mMD <- readODS::read_ods("~/proj5HT/den/serotonin/Data/Master_UI_data.ods", sheet = "Master_Macaque")
mMD <- mMD[!is.na(mMD$cell_name), ]
cells = unique(mMD$cell_name)



# cells already added to rookery
nnested_cells = sh5HT$nnested_cells

# macaque user input for experiments
output_cellsMD = MD[MD$cell_name %in% cells,]

#### Build ####
update_cells = cells[!cells %in% sh5HT$nnested_cells]

### Basic examples
i = "QF25.26.023.19.06.03"
i = "QF25.26.024.19.04.05"
i = "QF25.26.024.19.05.01"

## Way wash in
i = "QF25.26.024.19.03.03"
i = "QF25.26.024.19.04.03"
i = "QF25.26.024.19.08.03"
i = "QF25.26.024.19.08.03"
i = "QF25.26.016.11.01.02"
i =  "QF25.26.023.19.06.03"


cell = "QN26.26.004.20.05.01"
cell = "QN26.26.004.20.02.03"
cell = "QN26.26.004.20.03.03"
cell = "QN26.26.004.20.08.02"
cell = "QN25.26.017.20.06.04"
cell = "QM26.26.022.20.04.03"
cell = "QM26.26.022.20.06.02"
cell = "QN26.26.007.20.03.03"
cell = "QN22.26.011.1A.38.02"

cell = "H22.03.301.11.09.01.03"

#### Macaque cells
mMD <- readODS::read_ods("~/proj5HT/den/serotonin/Data/Master_UI_data.ods", sheet = "Master_Macaque")
mMD <- mMD[!is.na(mMD$cell_name), ]
cells = unique(mMD$cell_name)
output_cellsMD = MD[MD$cell_name %in% cells,]

#### Human cells
hMD <- readODS::read_ods("~/proj5HT/den/serotonin/Data/Master_UI_data.ods", sheet = "Human")
hMD <- hMD[!is.na(hMD$cell_name), ]
cells = unique(hMD$cell_name)
output_cellsMD = MD[MD$cell_name %in% cells,]

cell = cells[94]
cell = update_cells[9]

which(cell == update_cells)
which(cell == cells)

cell = cells[30]

cell = tail(update_cells, 3)[1]



for (cell in tail(update_cells, 3)) {
  message("\n===== ", cell, " =====")
    srt <- load5HT(cell, tag = "-srt.rds", rehydrate = FALSE)
    #file.remove(list.files(srt$rd, full.names = TRUE)[grep("standard_puff", list.files(srt$rd))])}

  tryCatch({


    list.files(srt$rd)

    if (is.null(srt)) {
      srt <- nnest5HT(cell)
    }
    if (!is.list(srt)) return(invisible(NULL))

    srt$dfs$asp = build_asp(srt, mMD = mMD, write_csv = TRUE)
    srt$dfs$rmp = build_rmp(srt, mMD = mMD, save_srt = FALSE, write_csv = TRUE)

    #plot_rmp_pulse_features(srt$dfs$rmp$by_protocol[[1]]$rmp_protocol_analysis$pulse)

    srt$dfs$blSpike = spikePuff_from_asp(srt, mMD = mMD, plot_it = TRUE)
    srt$dfs$dwin = DrugWashIn_from_asp(srt, mMD = mMD, plot_it = TRUE)



    if(!is.null(srt$dfs$asp)){
    srt = make_figs(srt, y_mode = "percent", x_lim_by_protocol = c(-15,50))
    #
    # plot_asp_single_protocol(
    #   srt$dfs$asp,
    #   protocol_key = srt$dfs$asp$protocol_keys[1],
    #   heat_dt = 1,
    #   ttl_shape = 6,
    #   y_mode = "log2_fc",
    #   #y_lim = y_lim_single,
    #   x_lim = c(-10,55),
    #   show_plot = TRUE,
    #   save_png = TRUE
    # )
    # plot_asp_single_protocol(
    #   srt$dfs$asp,
    #   protocol_key = srt$dfs$asp$protocol_keys[5],
    #   heat_dt = 1,
    #   ttl_shape = 6,
    #   y_mode = "log2_fc",
    #   #y_lim = y_lim_single,
    #   x_lim = c(-10,55),
    #   show_plot = TRUE,
    #   save_png = TRUE
    # )


    }

    #res <- run_asp_fig_qc_loop_perkey(srt)
    #srt <- res$srt
    #call <- res$call
    #calls <- dplyr::bind_rows(calls, tibble::tibble(cell = srt$cell, call = call))

    saveRDS(srt, file = srt$child_files$nnest)


  }, error = function(e) {
    # Save a debug bundle and stop (your preference)
    dbg_file <- file.path(tempdir(), paste0("spikePuff5HT_ERROR_", cell, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))
    saveRDS(list(cell = cell, error = e, traceback = sys.calls()), dbg_file)
    message("ERROR in cell ", cell, " — saved debug to: ", dbg_file)
    stop(e)  # stop immediately
  })

}

srt$figs$concat <- plot_asp_concatenated(
  srt$dfs$asp,
  gap_s = 2,
  heat_dt = 1,
  order_by = "wrapped_sweeps",
  y_mode = "log2_fc",
  y_lim = NULL,
  show_gap_lines = TRUE,
  save_png = TRUE,
  png_path = "figs/exp_concat_heat/",
  verbose = TRUE
)


cells_to_exclude = c(
  "Q21.26.010.11.13.03"
)


# Assemble complete dataframe ---------------------------------------------

comp5HT = compile5HT()
comp5HT = readRDS("data-raw/compiled5HT.rds")

#### archive ####

reticulate::source_python("py/plotting_pdf.py")

plot_nwb_sweeps(
  nwb_path = srt$files$nwb,
  sweeps   = list(40L),
  out_dir  = "figs/traces",
  fmt      = "pdf",
  y_unit_target = "mV",
  pub_void_axes = TRUE,
  scalebar_x = 10,
  scalebar_y = 10
)



#
# # ---- decision loop ----
# call <- NA_character_
#
# # ---- Figures plus call ----
# repeat {
#
#   srt <- make_figs(srt,
#                    y_lim_single = c(-1, NA),
#                    y_lim_concat = c(NA, NA))
#
#   call <- select.list(
#     choices = c("keep", "omit", "rerun_figures"),
#     title   = srt$cell
#   )
#
#   # user hit cancel / closed dialog
#   if (is.null(call) || identical(call, "")) {
#     message("No selection made; asking again.")
#     next
#   }
#
#   if (call %in% c("keep", "omit")) break
#
#   # otherwise call == "rerun_figures" -> loop continues
# }
#
# # record decision
# calls <- dplyr::bind_rows(
#   calls,
#   tibble::tibble(cell = srt$cell, call = call)
# )
#
# askYesNo("move to next")
# srtPlots = function(x){
#   tst = srt$dfs$spikeTTL$spike_puff_output
#   plot(tst$time, tst$percent_change, main = paste(srt$dfs$spikeTTL$selected_MD$cell_name, srt$dfs$spikeTTL$iMerged$predicted[!is.na(srt$dfs$spikeTTL$iMerged$predicted)]), xlab = "time", ylab = "% Change")
# }

# for(i in cells){
#   print(i)
#
#   #srt = nnest5HT(i)
#   srt = projHCT::loadHCT(i, tag = "-srt.rds")
#
#   if (is.null(srt)) {
#     srt = nnest5HT(i)
#     if(!file.exists(srt$files$hct_nnest))
#       srt = projHCT::nnestHCT(i)
#     srt = nnest5HT(i)
#   }
#   if(!is.list(srt)){
#     next
#   }
#
#   if("puff" %in% names(srt$dfs)){
#     srt$dfs$puff = NULL
#     saveRDS(srt, file = srt$files$nnest)
#   }
#
#   srt$dfs$spikeTTL = spikePuffDrug(srt, mMD = mstr)
#
#   srt$dfs$DrugWash = DrugWashIn(srt)
#
#   plotBasicPuff(srt)
#
#   plotBasicWash(srt)
#
#   saveRDS(srt, file = srt$files$nnest)
#
#
# }

#### Sodium availability ####
#NaAv = NaAvail5HT(srt, show_plot = TRUE, return_dfs = TRUE)
# NaAvail5HT(srt, show_plot = TRUE, return_dfs = FALSE)
#
#
# df_na = srt$dfs$NaAv$df_na
# df_spike= srt$dfs$NaAv$df_spike
#
# ggplot() +
#   # sodium availability background
#   geom_line(
#     data = df_na,
#     aes(x = time, y = a, color = cond),
#     linewidth = 0.6, alpha = 0.7
#   ) +
#   # spike %-change points (rescaled)
#   geom_point(
#     data = df_spike,
#     aes(x = time, y = percent_scaled),
#     shape = 21, fill = "black", color = "white",
#     size = 2, alpha = 0.9
#   ) +
#   # vertical marker for drug onset
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
#   annotate("text", x = 0, y = 1.02, label = "Drug onset", hjust = -0.1, size = 4) +
#   labs(
#     title = paste(
#       "Continuous Sodium Availability and Spike % Change —",
#       srt$dfs$spikeTTL$selected_MD$cell_name
#     ),
#     x = "Time (s, 0 = drug onset)",
#     y = "Sodium Availability a(t)"
#   ) +
#   scale_y_continuous(
#     name = "Sodium Availability a(t)",
#     sec.axis = sec_axis(
#       trans = ~ . * (max(df_spike$percent_change, na.rm = TRUE) -
#                        min(df_spike$percent_change, na.rm = TRUE)) +
#         min(df_spike$percent_change, na.rm = TRUE),
#       name = "% Change"
#     )
#   ) +
#   scale_color_manual(values = c("Baseline" = "salmon", "Drug" = "cyan3")) +
#   theme_minimal(base_size = 14) +
#   theme(
#     panel.grid.minor = element_blank(),
#     legend.position = "right",
#     plot.title = element_text(face = "bold")
#   )
#
# plot(gg)

