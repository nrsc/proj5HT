### Project looper
library(dplyr)
library(tidyr)
library(nphys)
library(projHCT)
library(proj5HT)

#### Setup ####
MD = projHCT::sheets$MD



#sheets = projHCT::sheets$files$nnest_files
sh5HT = sheets5HT()
#sh5HT = proj5HT::sh5HT

# cells already added to rookery
cells = sh5HT$cells
# macaque user input for experiments
sMdta = sh5HT$mUI
output_cellsMD = MD[MD$cell_name %in% cells,]
#cells = grep("QF25.26.022",cells, value = TRUE)

#lf = list.files("rookery")

#i = cells[10]
#### Build ####
# Want to build baseline appreciating function
#MD[grep(i, MD$cell_name),]

update_cells = sh5HT$mUI$cell_name[!sh5HT$mUI$cell_name %in% sh5HT$cells]

i = "QF25.26.023.19.06.03"
i = "QF25.26.016.11.01.02"


i = cells[1]
cells[21]
which(cells == "QN25.26.017.11.02.02")
which(cells == i)
which(update_cells == i)

i = update_cells[4]

for(i in update_cells){
  print(i)

  #srt = nnest5HT(i)
  srt = projHCT::loadHCT(i, tag = "-srt.rds")

  # if (!is.null(srt)) {
  # srt = nnest5HT(srt)
  # }

  if (is.null(srt)) {
      srt = nnest5HT(i)
    if(!file.exists(srt$files$hct_nnest))
      srt = projHCT::nnestHCT(i)
      srt = nnest5HT(i)
  }
  if(!is.list(srt)){
    next

  }

  srt$cell
  srt$dfs$puff = NULL

  srt$dfs$NaAv = NULL

  #file.exists(srt$files$nwb)



  #srt = spikePuff5HT(srt)
  #srt$dfs$spikeTTL$selected_MD

  #srt = ketWashIn(srt)

  #NaAv = NaAvail5HT(srt, show_plot = TRUE, return_dfs = TRUE)
  #srt$dfs$NaAv = NaAvail5HT(srt, show_plot = TRUE, return_dfs = TRUE)

  #if(file.exists(file.path(srt$rd, paste0(srt$cell, "-srtPuff.csv")))){
  #  srt$files$exp_dfs$srtPuff = file.path(srt$rd, paste0(srt$cell, "-srtPuff.csv"))
  #}
  #if(file.exists(file.path(srt$rd, paste0(srt$cell, "-srtKetWash.csv")))){
  #  srt$files$exp_dfs$ketWash = file.path(srt$rd, paste0(srt$cell, "-srtKetWash.csv"))
  #}

  #askYesNo()

  saveRDS(srt, file = srt$files$nnest)


}

# Assemble complete dataframe ---------------------------------------------

comp5HT = compile5HT()

srtPlots = function(x){

  tst = srt$dfs$spikeTTL$spike_puff_output
  plot(tst$time, tst$percent_change, main = paste(srt$dfs$spikeTTL$selected_MD$cell_name, srt$dfs$spikeTTL$iMerged$predicted[!is.na(srt$dfs$spikeTTL$iMerged$predicted)]), xlab = "time", ylab = "% Change")

  NaAvail5HT(srt, show_plot = TRUE, return_dfs = FALSE)


  df_na = srt$dfs$NaAv$df_na
  df_spike= srt$dfs$NaAv$df_spike

  ggplot() +
    # sodium availability background
    geom_line(
      data = df_na,
      aes(x = time, y = a, color = cond),
      linewidth = 0.6, alpha = 0.7
    ) +
    # spike %-change points (rescaled)
    geom_point(
      data = df_spike,
      aes(x = time, y = percent_scaled),
      shape = 21, fill = "black", color = "white",
      size = 2, alpha = 0.9
    ) +
    # vertical marker for drug onset
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    annotate("text", x = 0, y = 1.02, label = "Drug onset", hjust = -0.1, size = 4) +
    labs(
      title = paste(
        "Continuous Sodium Availability and Spike % Change —",
        srt$dfs$spikeTTL$selected_MD$cell_name
      ),
      x = "Time (s, 0 = drug onset)",
      y = "Sodium Availability a(t)"
    ) +
    scale_y_continuous(
      name = "Sodium Availability a(t)",
      sec.axis = sec_axis(
        trans = ~ . * (max(df_spike$percent_change, na.rm = TRUE) -
                         min(df_spike$percent_change, na.rm = TRUE)) +
          min(df_spike$percent_change, na.rm = TRUE),
        name = "% Change"
      )
    ) +
    scale_color_manual(values = c("Baseline" = "salmon", "Drug" = "cyan3")) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold")
    )

  plot(gg)



}
