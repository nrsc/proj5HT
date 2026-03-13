#### Figure 2 ####

df %>%
  dplyr::filter(
    assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET", "L5_PN"),
    bath == "Ket") %>%
  group_by(x_bin, cell_name, assigned_subclass) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  group_by(cell_name, assigned_subclass) %>%
  tidyr::complete(
    x_bin = seq(min(x_bin), max(x_bin), by = 5),
    fill = list(y = 0)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = x_bin, y = y, colour = cell_name)) +
  geom_line() +
  facet_wrap(~assigned_subclass) +
  ylim(0, 200) +
  xlim(-20, 100) +
  ylab("Percent change") +
  xlab("Time (s)") +
  theme_minimal() +
  theme(legend.position = "none")

