plotBasicWash = function(x){

  if(!"DrugWash" %in% names(x$dfs)){
    message("No drug washin to plot")
    return(NULL)
  }

  df = x$dfs$spikeTTL$spike_puff_output
  dw = x$dfs$DrugWash$spike_puff_output
  df = rbind(df,dw)

  title_string = paste(x$md$cell_name, unique(df$assigned_subclass), unique(df$puff), "Wash",  unique(dw$bath), sep = "-")


  df_binned <- df %>%
    mutate(time_sec = floor(time)) %>%
    group_by(time_sec, bath) %>%
    summarise(
      percent_change = mean(percent_change, na.rm = TRUE),
      instRate = mean(instRate, na.rm = TRUE),
      .groups = "drop"
    )

  # 3. Plot with dual axis
  gg = ggplot(df_binned, aes(x = time_sec)) +
    # Left axis: % change
    geom_line(data = df_binned, aes(time_sec, percent_change, colour = bath), linewidth = 1) +
    geom_point(aes(y = percent_change, colour = bath)) +
    scale_colour_manual(values = c("blue", "red")) +
    xlim(-20,60) +
    labs(
      x = "Time (sec)",
      title = title_string
    ) +
    theme_bw()


  plot(gg)

  ggsave(file.path(srt$rd, paste0(title_string, ".png")), plot = gg)
  ggsave(file.path("figs/washIn", paste0(title_string, ".png")), plot = gg)

  gg2 = ggplot(df_binned, aes(x = time_sec)) +
    # Left axis: % change
    geom_line(data = df_binned, aes(time_sec, instRate, colour = bath), linewidth = 1) +
    geom_point(aes(y = instRate, colour = bath)) +
    scale_colour_manual(values = c("blue", "red")) +
    xlim(-20,60) +
    labs(
      x = "Time (sec)",
      title = title_string
    ) +
    theme_bw()


  plot(gg2)

  ggsave(file.path(srt$rd, paste0(title_string,"-instRate", ".png")), plot = gg2)
  ggsave(file.path("figs/washIn", paste0(title_string,"-instRate", ".png")), plot = gg2)


}
