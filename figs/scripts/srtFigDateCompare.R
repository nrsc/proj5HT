dt = as.Date("2025-11-16")
Before = paste("Before", dt)
After = paste("After", dt)

df = df %>% mutate(
  date_clean = as.Date(stringr::str_remove(Date, "\\s+PZT")),
  date_group = if_else(
    date_clean < dt,
    Before,
    After
  ),
  date_group = factor(date_group, levels = c(Before, After))
)


# All in one --------------------------------------------------------------
df %>%
  dplyr::filter(
    assigned_subclass %in% c("L5_ET", "L5_IT"),
    #assigned_subclass %in% c("L23_IT"),
    bath == "none",
    expCon == "Standard_Puff",
    puff %in% c("5HT[100uM]", "5HT[200uM]")
  ) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species, date_group) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(aes(
    x = x_bin,
    y = y,
    group = cell_name,
    colour = assigned_subclass
  )) +
  geom_line(alpha = 0.5) +
  facet_grid(
    rows = vars(Species),
    cols = vars(assigned_subclass, date_group)
  ) +
  ylim(0, 210) +
  xlim(-50, 100) +
  ylab("Percent Change") +
  xlab("Time (s)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 9)
  )

df %>%
  dplyr::filter(
    assigned_subclass %in% c("L23_IT"),
    bath == "none",
    expCon == "Standard_Puff",
    puff %in% c("5HT[100uM]", "5HT[200uM]")
  ) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species, date_group) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(aes(
    x = x_bin,
    y = y,
    group = cell_name,
    colour = assigned_subclass
  )) +
  geom_line(alpha = 0.5) +
  facet_grid(
    rows = vars(Species),
    cols = vars(assigned_subclass, date_group)
  ) +
  ylim(0, 320) +
  xlim(-50, 100) +
  ylab("Percent Change") +
  xlab("Time (s)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 9)
  )

df %>%
  dplyr::filter(
    assigned_subclass %in% c("L5_ET", "L5_IT"),
    #assigned_subclass %in% c("L23_IT"),
    bath == "none",
    expCon == "Standard_Puff",
    puff %in% c("5HT[100uM]", "5HT[200uM]")
  ) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species, date_group) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(aes(
    x = x_bin,
    y = y,
    group = cell_name,
    colour = assigned_subclass
  )) +
  geom_line(alpha = 0.5) +
  facet_grid(
    #rows = vars(Species),
    cols = vars(assigned_subclass, date_group)
  ) +
  ylim(0, 250) +
  xlim(-50, 100) +
  ylab("Percent Change") +
  xlab("Time (s)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 9)
  )

plot(gg1)
