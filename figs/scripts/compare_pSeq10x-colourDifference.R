df_diag <- df_diag %>%
  dplyr::mutate(delta_log1p_cpm = log1p(cpm_brain) - log1p(cpm_cMTG))

lim <- max(abs(df_diag$delta_log1p_cpm), na.rm = TRUE)  # 2.3100432

ggplot2::ggplot(df_diag, ggplot2::aes(gene, subclass)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_brain, color = delta_log1p_cpm)) +
  ggplot2::geom_point(ggplot2::aes(size = pct_cMTG), shape = 1, stroke = 1) +
  ggplot2::scale_color_gradient2(
    limits = c(-lim, lim),
    midpoint = 0,
    low = "yellow",
    mid = "grey15",
    high = "chartreuse3",
    oob = scales::squish
  ) +
  ggplot2::labs(color = "Δ log1p(CPM)\n(pSeq - 10x)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
