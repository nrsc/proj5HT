plotBasicPuff = function(x){

  if(!"spikeTTL" %in% names(x$dfs)){
    message("SpikeTTL to plot")
    return(NULL)
  }

  df = x$dfs$spikeTTL$spike_puff_output

title_string = paste(x$md$cell_name, unique(df$assigned_subclass), unique(df$puff), "BasicPuff", sep = "-")

# 1. Bin into 1-second bins
df_binned <- df %>%
  mutate(time_sec = floor(time)) %>%
  group_by(time_sec) %>%
  summarise(
    percent_change = mean(percent_change, na.rm = TRUE),
    instRate = mean(instRate, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Compute scaling between percent_change and instRate
pc_range  <- range(df_binned$percent_change, na.rm = TRUE)
ir_range  <- range(df_binned$instRate, na.rm = TRUE)

# scale factor to map instRate → percent_change-axis
scale_factor <- diff(pc_range) / diff(ir_range)

df2 <- df %>% mutate(time_sec = floor(time))

# 3. Plot with dual axis
gg = ggplot(df_binned, aes(x = time_sec)) +
  # Left axis: % change
  geom_line(data = df_binned, aes(time_sec, percent_change), color = "red", linewidth = 1) +
  geom_point(aes(y = percent_change)) +
  geom_point(data = df2, aes(time, percent_change), alpha = 0.3) +

  # Right axis: instRate (scaled)
  geom_line(aes(y = instRate * scale_factor),
            color = "blue", linewidth = 1) +

  scale_y_continuous(
    name = "% change",
    sec.axis = sec_axis(~ . / scale_factor, name = "Inst. Rate (Hz)")
  ) +
  labs(
    x = "Time (sec)",
    title = paste(df$cell_name, df$assigned_subclass, df$puff, sep = " - ")
  ) +
  theme_bw()

plot(gg)

ggsave(file.path(srt$rd, paste0(title_string, ".png")), plot = gg)
ggsave(file.path("figs/basicPuff", paste0(title_string, ".png")), plot = gg)


}
