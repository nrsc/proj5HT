figures5HT = function(){

  #comp5HT = compile5HT()
  df = data.frame(comp5HT$srtPuff)#, comp5HT$srtKetWash)

  unique(df$expCon)
  unique(df$puff)
  unique(df$bath)
  df$Species



  df %>% subset(., assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET", "L5_PN") & bath == "none" & expCon == "Standard_Puff" & puff == "5HT[100uM]") %>%
    #lfout %>% subset(., assigned_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET", "L5_PN") & ket == FALSE) %>%
    #lfout %>% subset(., !ket) %>%
    mutate(x_bin = floor(time / 5) * 5) %>%
    group_by(x_bin, cell_name, assigned_subclass, Species) %>%
    summarise(y = mean(percent_change), .groups = "drop") %>%
    ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
    #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
    geom_line() +
    geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
    facet_wrap(~assigned_subclass) +
    facet_grid(~Species) +
    #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
    ylim(0, 500) +
    ylab("Time (s)") +
    xlim(-50, 100) +
    theme_minimal() +
    theme(legend.position = "none")


  #
  df %>% subset(., assigned_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET", "L5_PN") & ket == FALSE & expCon == "Standard_Puff") %>%
    #lfout %>% subset(., assigned_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET", "L5_PN") & ket == FALSE) %>%
    #lfout %>% subset(., !ket) %>%
    mutate(x_bin = floor(time / 5) * 5) %>%
    group_by(x_bin, cell_name, assigned_subclass) %>%
    summarise(y = mean(percent_change), .groups = "drop") %>%
    ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
    #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
    geom_line() +
    geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
    facet_wrap(~assigned_subclass) +
    #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
    ylim(0, 500) +
    ylab("Time (s)") +
    xlim(-50, 100) +
    theme_minimal() +
    theme(legend.position = "none")
  #

  df %>% subset(., assigned_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET", "L5_PN") & expCon == "Ket" | expCon == "Ket-wash-in") %>%
    #lfout %>% subset(., !ket) %>%
    mutate(x_bin = floor(time / 5) * 5) %>%
    #mutate(x_bin = cut(time, breaks = seq(0, max(time) + 5, by = 5), right = FALSE)) %>%
    group_by(x_bin, cell_name, assigned_subclass) %>%
    summarise(y = mean(percent_change), .groups = "drop") %>%
    ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
    #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
    geom_line() +
    geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
    facet_wrap(~assigned_subclass) +
    ggtitle("Ketanserin") +
    #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
    ylim(0, 200) +
    ylab("Time (s)") +
    xlim(-20, 100) +
    theme_minimal() +
    theme(legend.position = "none")
  #
  df %>% subset(., assigned_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET", "L5_PN") & expCon == "Control") %>%
    #lfout %>% subset(., !ket) %>%
    mutate(x_bin = floor(time / 5) * 5) %>%
    #mutate(x_bin = cut(time, breaks = seq(0, max(time) + 5, by = 5), right = FALSE)) %>%
    group_by(x_bin, cell_name, assigned_subclass) %>%
    summarise(y = mean(percent_change), .groups = "drop") %>%
    ggplot(., aes(x = x_bin, y = y, colour = cell_name)) +
    #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
    ggtitle("control experiments") +
    geom_line() +
    geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
    facet_wrap(~assigned_subclass) +
    #facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
    ylim(0, 200) +
    ylab("Time (s)") +
    xlim(-50, 50) +
    theme_minimal() +
    theme(legend.position = "none")




# # Predicted Subclass plots ------------------------------------------------
# lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
#   mutate(x_bin = floor(t / 5) * 5) %>%
#   group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
#   summarise(y = mean(percent_change), .groups = "drop") %>%
#   ggplot(., aes(x = x_bin, y = y, colour = cell)) +
#   #ggplot(., aes(x = t, y = rate, colour = cell)) +
#   geom_line() +
#   geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
#   facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
#   ggtitle("% change grouped by predicted_subclass") +
#   ylim(0, 300) +
#   ylab("Time (s)") +
#   #xlim(-50, 120) +
#   xlim(-25, 60) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
#
# lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
#   mutate(x_bin = floor(t / 5) * 5) %>%
#   group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
#   summarise(y = mean(rate), .groups = "drop") %>%
#   ggplot(., aes(x = x_bin, y = y, colour = cell)) +
#   #ggplot(., aes(x = t, y = rate, colour = cell)) +
#   geom_line() +
#   ggtitle("Firing rate grouped by predicted_subclass") +
#   geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
#   facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
#   #ylim(0, 300) +
#   ylab("Time (s)") +
#   xlim(-50, 120) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
#
#
# # Subclass_Tree label plots -----------------------------------------------
#
# #### Human vs Macacque subclass_Tree % Change
# lfout %>% subset(., subclass_Tree %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
#   mutate(x_bin = floor(t / 5) * 5) %>%
#   group_by(x_bin, protocol, cell, subclass_Tree, Species) %>%
#   summarise(y = mean(percent_change), .groups = "drop") %>%
#   ggplot(., aes(x = x_bin, y = y, colour = cell)) +
#   #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
#   geom_line() +
#   geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
#   facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
#   ylim(0, 300) +
#   ylab("Time (s)") +
#   xlim(-50, 120) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
#
#
# #### Human vs Macacque subclass_Tree rate
# lfout %>% subset(., subclass_Tree %in% c("L2/3_IT", "L5_IT", "L5_ET")) %>%
#   mutate(x_bin = floor(t / 5) * 5) %>%
#   group_by(x_bin, protocol, cell, subclass_Tree, Species) %>%
#   summarise(y = mean(rate), .groups = "drop") %>%
#   ggplot(., aes(x = x_bin, y = y, colour = cell)) +
#   #ggplot(., aes(x = t, y = rate, colour = cell)) +
#   geom_line() +
#   geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
#   facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
#   #ylim(0, 300) +
#   ylab("Time (s)") +
#   xlim(-50, 120) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
#
#
# # archive -----------------------------------------------------------------
# lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
#   mutate(x_bin = floor(t / 5) * 5) %>%
#   group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
#   summarise(y = mean(percent_change), .groups = "drop") %>%
#   ggplot(., aes(x = x_bin, y = y, colour = cell)) +
#   geom_line() +
#   #facet_wrap( ~ predicted_subclass, nrow = 3, scales = "free") +
#   facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
#   ylim(0, 300) +
#   ylab("Time (s)") +
#   xlim(-25, 120) +
#   theme_minimal() +
#   theme(legend.position = "none")

}
