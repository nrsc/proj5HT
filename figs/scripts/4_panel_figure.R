plot_response_types_4panel <- function(df_prepped,
                                       per_cell,
                                       xlim = c(-5, 20),
                                       k_smooth = 11,
                                       alpha = 0.35) {

  keep_tbl <- per_cell %>%
    dplyr::mutate(
      resp_class = dplyr::case_when(
        response_type == "excitation" ~ "Excitation",
        response_type == "inhibition" ~ "Inhibition",
        response_type == "biphasic" & biphasic_order == "exc_then_inh" ~ "Biphasic exc→inh",
        response_type == "biphasic" & biphasic_order == "inh_then_exc" ~ "Biphasic inh→exc",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(resp_class)) %>%
    dplyr::select(cell_name, resp_class) %>%
    dplyr::distinct() %>%
    dplyr::mutate(
      resp_class = factor(resp_class,
                          levels = c("Excitation","Inhibition","Biphasic exc→inh","Biphasic inh→exc"))
    )

  d <- df_prepped %>%
    dplyr::inner_join(keep_tbl, by = "cell_name") %>%
    dplyr::mutate(time = as.numeric(time)) %>%
    dplyr::filter(time >= xlim[1], time <= xlim[2]) %>%
    dplyr::group_by(resp_class, cell_name) %>%
    dplyr::arrange(time, .by_group = TRUE) %>%
    dplyr::mutate(
      pct_smooth = zoo::rollmean(pct_of_baseline, k = k_smooth,
                                 fill = NA, align = "center")
    ) %>%
    dplyr::ungroup()

  cols <- c(
    "Excitation" = "#d62728",
    "Inhibition" = "#1f77b4",
    "Biphasic exc→inh" = "purple",
    "Biphasic inh→exc" = "darkgreen"
  )

  ggplot(d, aes(x = time, y = pct_smooth, group = cell_name)) +
    facet_wrap(~resp_class, ncol = 2) +
    geom_hline(yintercept = 100, linetype = 2, linewidth = 0.5, color = "grey40") +
    geom_line(aes(color = resp_class), alpha = alpha, linewidth = 0.7, na.rm = TRUE) +
    scale_color_manual(values = cols, guide = "none") +
    labs(
      title = "Response types (4-panel view)",
      subtitle = "Smoothed %baseline traces; each panel is one response class",
      x = "Time (s)",
      y = "% baseline"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())
}

# Example:
p4 <- plot_response_types_4panel(df_prepped, per_cell2, xlim = c(-5, 50), k_smooth = 3)
print(p4)
