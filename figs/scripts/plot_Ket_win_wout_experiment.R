srt <- load5HT("QN26.26.004.20.03.03", tag = "-srt.rds")

res1 <- plot_asp_single_protocol(
  srt$dfs$asp,
  protocol_key = names(srt$dfs$asp$by_protocol)[1],
  heat_dt = 0.5,
  ttl_shape = 6,
  y_mode = "log2_fc",
  y_lim = c(-.5, 1)
)
res5 <- plot_asp_single_protocol(
  srt$dfs$asp,
  protocol_key = names(srt$dfs$asp$by_protocol)[5],
  heat_dt = 0.5,
  ttl_shape = 6,
  y_mode = "log2_fc",
  y_lim = c(-.5, 1)
)
res8 <- plot_asp_single_protocol(
  srt$dfs$asp,
  protocol_key = names(srt$dfs$asp$by_protocol)[8],
  heat_dt = 0.5,
  ttl_shape = 6,
  y_mode = "log2_fc",
  y_lim = c(-.5, 1)
)

r3 <- (res1$plot | res5$plot | res8$plot) +
  patchwork::plot_layout(widths = c(1, 1, 1)) &
  theme(
    plot.subtitle = element_blank(),   # ensure no subtitles anywhere
    plot.title = element_text(size = 11, face = "bold"),
    plot.title.position = "plot"
  )

