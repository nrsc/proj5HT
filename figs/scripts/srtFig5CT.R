dfs$df_5ct %>% #subset(., assigned_subclass %in% c("L5_IT", "L5_ET")) %>%
  #mutate(x_bin = floor(time / 2) * 2) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, group = cell_name, colour = assigned_subclass)) +
  #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
  geom_line() +
  #geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_wrap(~assigned_subclass) +
  #facet_grid(~Species) +
  #facet_grid(rows = vars(Species), cols = vars(assigned_subclass)) +
  ylim(0, 200) +
  ggtitle("5CT response across macaque pyramids") +
  xlab("Time (s)") +
  ylab("Percent Change") +
  xlim(-10, 50) +
  theme_minimal() +
  theme(legend.position = "none")
