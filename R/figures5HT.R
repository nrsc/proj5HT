figures5HT = function(){

  #### Import dataframes ####
  source("~/proj5HT/R/figures_df_build.R")
  #file.edit("~/proj5HT/R/figures_df_build.R")
  dfs = figure_df_build()

  unique(dfs$dfH$x_bin)
  range(dfs$dfH$percent_change)

  #### Human data ####
  #file.edit("~/proj5HT/figs/scripts/srtHumanData.R")
  source("~/proj5HT/figs/scripts/srtHumanData.R", echo = TRUE)


  #### Figure comparing dates of L5 and L2/3 cells by date ####
  #file.edit("figs/scripts/srtFigDateCompare.R")
  file.edit("figs/scripts/srt_ketContam.R")
  source("figs/scripts/srtFigDateCompare.R", echo = TRUE)





  #### Ket figures for Macaque ####
  #file.edit("~/proj5HT/figs/scripts/srtFigKet.R")
  source("~/proj5HT/figs/scripts/srtFigKet.R", echo = TRUE)


  #### Ket figures for Macaque ####
  #file.edit("~/proj5HT/figs/scripts/srtFig5CT.R")
  source("~/proj5HT/figs/scripts/srtFig5CT.R", echo = TRUE)




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



  df_ket %>% subset(., assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET") & bath == "Ket" & expCon == "Standard_Puff" &  puff %in% c("5HT[100uM]", "5HT[200uM]")) %>%
    #mutate(x_bin = floor(time / 5) * 5) %>%
    group_by(x_bin, cell_name, assigned_subclass) %>%
    summarise(y = mean(percent_change), .groups = "drop") %>%
    ggplot(., aes(x = x_bin, y = y, group = cell_name, colour = assigned_subclass)) +
    #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
    geom_line() +
    #geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
    facet_wrap(~assigned_subclass) +
    #facet_grid(~Species) +
    ##facet_grid(rows = vars(Species), cols = vars(assigned_subclass)) +
    #ylim(0, 500) +
    ylab("Time (s)") +
    xlab("Percent Change") +
    xlim(-10, 40) +
    theme_minimal() +
    theme(legend.position = "none")






  dfs$df %>% subset(., assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET", "L5_PN") & expCon == "Control") %>%
    #lfout %>% subset(., !ket) %>%
    #mutate(x_bin = floor(time / 5) * 5) %>%
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
    ylim(50, 150) +
    ylab("Time (s)") +
    xlim(-50, 50) +
    theme_minimal() +
    theme(legend.position = "none")

  file.edit("~/proj5HT/figs/scripts/srtFig5CTpuff.R")

  dfW = data.frame(comp5HT$srtDrugWash)

}
