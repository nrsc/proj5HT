#### Figure 2 ####
#comp5HT = compile5HT()
df = data.frame(comp5HT$srtPuff)
# This figure will show the basic 5HT puff protocol across different cell types
#### All cells types and species
df %>% subset(., assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET", "L5_PN") & bath == "none" & expCon == "Standard_Puff" & puff == "5HT[100uM]") %>%
  mutate(x_bin = floor(time / 5) * 5) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
  #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
  geom_line() +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  #facet_wrap(~assigned_subclass) +
  facet_grid(rows = vars(Species), cols = vars(assigned_subclass)) +
  #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
  ylim(0, 500) +
  ylab("Time (s)") +
  xlim(-20, 60) +
  theme_minimal() +
  theme(legend.position = "none")

df %>% subset(., assigned_subclass %in% c("L5_ET") & bath == "none" & expCon == "Standard_Puff" & puff == "5HT[100uM]") %>%
  mutate(x_bin = floor(time / 5) * 5) %>%
  group_by(x_bin, cell_name, assigned_subclass, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
  #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
  geom_line() +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  #facet_wrap(~assigned_subclass) +
  facet_grid(rows = vars(Species), cols = vars(assigned_subclass)) +
  #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
  ylim(0, 500) +
  ylab("Time (s)") +
  xlim(-20, 60) +
  theme_minimal() +
  theme(legend.position = "none")
