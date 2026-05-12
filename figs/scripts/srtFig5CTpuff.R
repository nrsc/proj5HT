dfs$df_5ct %>% subset(., bath == "none") %>%
dfs$df_5ct %>% #subset(., bath == "none") %>%
  group_by(x_bin, cell_name, assigned_subclass) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
  geom_line() +
  #geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_wrap(~assigned_subclass) +
  #facet_grid(rows = vars(Species), cols = vars(assigned_subclass)) +
  #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
  ylim(0, 200) +
  ylab("Time (s)") +
  xlim(-10, 50) +
  theme_minimal()# +
  theme(legend.position = "none")
