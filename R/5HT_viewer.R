srt_viewer = function(x){

  srt = projHCT::loadHCT(x, tag = "-srt.rds")

  plot(
    tst$time,
    tst$percent_change,
    main = paste(
      srt$dfs$spikeTTL$selected_MD$cell_name,
      srt$dfs$spikeTTL$iMerged$predicted[!is.na(srt$dfs$spikeTTL$iMerged$predicted)]
    ),
    xlab = "time",
    ylab = "% Change"
  )


}
