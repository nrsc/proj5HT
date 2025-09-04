library(dplyr)
library(tidyr)

in_range <- function(x, lower, upper, inclusive = TRUE) {
  if (inclusive) {
    x >= lower & x <= upper
  } else {
    x > lower & x < upper
  }
}

### Script to present output from Brians python analysis
MD = projHCT::sheets$MD
cells = list.files("den/serotonin/output/")
cell_json = list.files("den/serotonin/output/", recursive = TRUE, full.names = TRUE, pattern = ".json")
merged_metadata = read.csv("den/serotonin/Data/merged_metadata.csv")
output_cellsMD = MD[MD$cell_name %in% cells,]
#length(output_cellsMD$cell_name)

#which
unique(merged_metadata$cell_name) %in% cells
output_cellsMD$cell_name %in% cells
length(output_cellsMD$cell_name %in% cells)

i = cells[12]

i = "QN25.26.003.12.04.02" #Check appropriate labeling of

i = "H21.03.309.11.06.01.01"
i = "H22.03.301.11.09.01.03"


which(cells == i)
cells = cells[-40]
cells[40]

# Want to build baseline appreciating function
#[30:length(cells)]
for(i in cells){
  print(i)

  srt = nnest5HT(i)

  srt$files$nnest

  iMD = MD[grep(i, MD$cell_name),]


  iMerged = merged_metadata[grep(i, merged_metadata$cell_name),]
  #srt = loadHCT(i)
  #smry = srt$exp$smry

  jfile = cell_json[grep(i, cell_json)]
  jlist = jsonlite::fromJSON(jfile)
  # If length of spikes for analysis is greater that 1
  if(length(jlist) > 1){
    spike_puff_sweeps = iMerged[which(iMerged$stimulus_description == "spike_puff"), "wrapped_sweeps"]
    spike_puff_sweeps = gsub("\\[|\\]", "", spike_puff_sweeps)
    spike_puff_sweeps = as.numeric(strsplit(spike_puff_sweeps, ",")[[1]])

    tst = lapply(1:length(jlist), function(x){
      jlist[[x]]$sweeps_analyzed %in% spike_puff_sweeps
    })

    n = which(sapply(tst, any))

  }else{
    n = 1
  }

  fp = data.frame(epoch_t = jlist[[n]]$spike_times, rate = jlist[[n]]$Instantaneous_frequency)


# Assessment of stimulus trace ----------------------------------------------------
  bin_width = 1

  ## Fix to any issues where analysis cuts off because no spiking happens after puff.
  maxBin = max(fp$epoch_t)
  if(maxBin < bin_width)
    maxBin = 50

  lvls = levels(cut(fp$epoch_t, breaks = seq(0, maxBin, by = bin_width), right = FALSE))
  df0 = data.frame(bin = lvls)

 out1 = fp %>% dplyr::mutate(bin = cut(fp$epoch_t, breaks = seq(0, maxBin, by = bin_width), right = FALSE, na.rm = FALSE)) %>%
    group_by(bin) %>%
    summarise(count = n(),
              rate = count / bin_width,
              .groups = "keep") %>%
    as.data.frame(.)

  out1 = merge(out1, df0, by = "bin", all = TRUE)

  if(nrow(out1) > 1 && is.na(out1$bin[nrow(out1)]))
    out1 = out1[-nrow(out1),]

  out1 <- out1 %>% mutate(across(everything(), ~tidyr::replace_na(.x, 0)))

  out1$t = seq(0, length.out = length(out1$bin), by = bin_width)

  pCh = mean(as.numeric(out1[which(out1$t < 5), "rate"]))

  out1$percent_change = out1$rate/pCh*100
  out1$protocol = "Stimulus"

  if(max(out1$percent_change) <= 200){
    ymax = 200
  }else{
    ymax = max(out1$percent_change)
  }

  #$gg = ggplot(df_binned, aes(x = t, y = percent_change)) +
  # gg = out1 %>% mutate(x_bin = floor(t / 5) * 5) %>%
  #   group_by(x_bin) %>%
  #   summarise(y = mean(percent_change), .groups = "drop") %>%
  #   ggplot(aes(x_bin, y)) +
  #     geom_point() +
  #     ggtitle(paste(iMD$cell_name, iMD$predicted_subclass, sep = "--")) +
  #     theme_minimal()


  # gg = ggplot(df_binned, aes(x_bin, y)) +
  #   geom_point() +
  #   ylim(0, ymax) +
  #   ggtitle(paste(iMD$cell_name, iMD$predicted_subclass, sep = "--")) +
  #   theme_minimal()

  #plot(gg)

  #ggsave(paste0("figs/pch5HT/", iMD$cell_name, "-stim5HT_pchPlot.png"), plot = gg)



# Baseline section --------------------------------------------------------

  if("baseline_spike_puff" %in% iMerged$stimulus_description) {
    bl_puff_sweeps = iMerged[which(iMerged$stimulus_description == "baseline_spike_puff"), "wrapped_sweeps"]
    bl_puff_sweeps = gsub("\\[|\\]", "", bl_puff_sweeps)
    bl_puff_sweeps = as.numeric(strsplit(bl_puff_sweeps, ",")[[1]])

    sa = jlist[[n]]$sweeps_analyzed[1]
    bin_width = 1

    if (in_range(bl_puff_sweeps, sa - 5, sa, inclusive = FALSE)) {
      sN = bl_puff_sweeps + 1
      nwb = extractNWB(iMD$nwb_file, sweeps = sN)
      #nwb$sampling_rate[sN]
      #plot(nwb$sweeps, type = "l")
      fpbl = APanalysis(nwb$sweeps, rate = nwb$sampling_rate[sN])
      fpbl = as.data.frame(fpbl)

      out2 = fpbl %>% dplyr::mutate(bin = cut(
        -fpbl$epoch_t,
        breaks = seq(floor(-max(fpbl$epoch_t)), 0, by = bin_width),
        right = TRUE,
        na.rm = FALSE
      )) %>%
        group_by(bin) %>%
        summarise(count = n(), rate = count / bin_width) %>%
        as.data.frame(.)

      out2$t = seq(-1, by = -bin_width, length.out = nrow(out2))
      out2$percent_change = out2$rate / pCh * 100

      out2$protocol = "Baseline"

      out2 %>% mutate(x_bin = floor(t / 5) * 5) %>%
        group_by(x_bin) %>%
        summarise(y = mean(percent_change), .groups = "drop") %>%
        ggplot(aes(x_bin, y)) +
        geom_point() +
        ylim(0, ymax) +
        ggtitle(paste(iMD$cell_name, iMD$predicted_subclass, sep = "--")) +
        theme_minimal()

      df = rbind(out2, out1)

      # tst0 = df %>% mutate(x_bin = floor(t / 5) * 5) %>%
      #   group_by(x_bin, protocol) %>%
      #   summarise(y = mean(percent_change), .groups = "drop")# %>%

      gg2 = df %>% mutate(x_bin = floor(t / 5) * 5) %>%
        group_by(x_bin, protocol) %>%
        summarise(y = mean(percent_change), .groups = "drop") %>%
        ggplot(aes(x_bin, y, colour = protocol)) +
        geom_point() +
        ylab("% Change") +
        xlab("Time(s)") +
        ylim(0, ymax) +
        ggtitle(paste(iMD$cell_name, iMD$predicted_subclass, sep = "--")) +
        theme_minimal()

      plot(gg2)

      ggsave(paste0("figs/pch5HT/", iMD$cell_name, "-blstim5HT_pchPlot.png"),
             plot = gg2)

    }
  }

# Build dataframe
df$cell_name = iMD$cell_name
df$predicted_subclass = iMD$predicted_subclass
df$subclass_Corr = iMD$subclass_Corr
df$subclass_Tree = iMD$subclass_Tree
df$Species = iMD$Species
df$Cortical_area = iMD$Cortical_area

  if(!file.exists(file.path("rookery", iMD$cell_name))){
    dir.create(file.path("rookery", iMD$cell_name))
  }

  write.csv(df, file.path("rookery", iMD$cell_name, paste0(iMD$cell_name, "-blStim-50sPuff.csv")))

}



# Assemble complete dataframe ---------------------------------------------

lf = list.files("rookery", pattern = "blStim-50sPuff.csv", recursive = TRUE, full.names = TRUE)

lfout = lapply(lf, function(l) {
  print(l)
  #r = read.csv(l)
  r = readr::read_csv(l, show_col_types = FALSE)
  if(nrow(r)==0){
    return(NA)
  }
  return(r)
}) %>% .[!is.na(.)] %>% bind_rows(.)

lfout$...1 = NULL
lfout$subclass_Corr = NULL
lfout$cell_name = lfout$cell

lfout =lfout %>% mutate(Species = if_else(Species == "Human", "Human", "Macaque"))
lfout =lfout %>% mutate(class_origin = if_else(subclass_Tree == "none", "classifier", "omics"))

lfout <- lfout %>%
  left_join(
    output_cellsMD %>% select(cell_name, Cortical_area),
    by = "cell_name"
  )
unique(lfout$Cortical_area)
lfout =lfout %>% mutate(Cortical_area = if_else(Cortical_area == "MTG", "TCx", Cortical_area))
unique(lfout$Cortical_area)
unique(lfout$cell_name) %in% unique(output_cellsMD$cell_name)

# Plots -------------------------------------------------------------------
## All cell types
# ggplot(lfout, aes(x = t, y = percent_change, colour = cell)) +
#   geom_line() +
#   facet_wrap(~predicted_subclass) +
#   ylim(0, max(lfout$percent_change)) +
#   xlim(0, 50) +
#   theme_minimal() +
#   theme(legend.position = "none")
#
# ggplot(subset(lfout, predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET")), aes(x = t, y = percent_change, colour = cell)) +
#   geom_line() +
#   facet_wrap(~predicted_subclass, nrow = 3, scales = "free") +
#   ggtitle("Pooled Species") +
#   xlim(-50, 120) +
#   theme_minimal() +
#   theme(legend.position = "none")


### Human vs Macaque % change
# ggplot(subset(lfout, predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Species %in% c("Human", "Macaque")), aes(x = t, y = percent_change, colour = cell)) +
#   geom_line() +
#   #facet_wrap(~predicted_subclass, nrow = 3, scales = "free") +
#   facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
#   ggtitle("Human vs Macaque") +
#   #ylim(0, 375) +
#   xlim(0, 50) +
#   theme_minimal() +
#   theme(legend.position = "none")


# Predicted Subclass plots ------------------------------------------------
lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
  mutate(x_bin = floor(t / 5) * 5) %>%
  group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell)) +
  #ggplot(., aes(x = t, y = rate, colour = cell)) +
  geom_line() +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
  ggtitle("% change grouped by predicted_subclass") +
  ylim(0, 300) +
  ylab("Time (s)") +
  #xlim(-50, 120) +
  xlim(-25, 60) +
  theme_minimal() +
  theme(legend.position = "none")


lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
  mutate(x_bin = floor(t / 5) * 5) %>%
  group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
  summarise(y = mean(rate), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell)) +
  #ggplot(., aes(x = t, y = rate, colour = cell)) +
  geom_line() +
  ggtitle("Firing rate grouped by predicted_subclass") +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
  #ylim(0, 300) +
  ylab("Time (s)") +
  xlim(-50, 120) +
  theme_minimal() +
  theme(legend.position = "none")



# Subclass_Tree label plots -----------------------------------------------

#### Human vs Macacque subclass_Tree % Change
lfout %>% subset(., subclass_Tree %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
  mutate(x_bin = floor(t / 5) * 5) %>%
  group_by(x_bin, protocol, cell, subclass_Tree, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell)) +
  #ggplot(., aes(x = t, y = percent_change, colour = cell)) +
  geom_line() +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
  ylim(0, 300) +
  ylab("Time (s)") +
  xlim(-50, 120) +
  theme_minimal() +
  theme(legend.position = "none")



#### Human vs Macacque subclass_Tree rate
lfout %>% subset(., subclass_Tree %in% c("L2/3_IT", "L5_IT", "L5_ET")) %>%
  mutate(x_bin = floor(t / 5) * 5) %>%
  group_by(x_bin, protocol, cell, subclass_Tree, Species) %>%
  summarise(y = mean(rate), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell)) +
  #ggplot(., aes(x = t, y = rate, colour = cell)) +
  geom_line() +
  geom_smooth(method = "loess", span = 0.2, se = FALSE, color = "red", size = 1.2) +
  facet_grid(cols = vars(Species), rows = vars(subclass_Tree)) +
  #ylim(0, 300) +
  ylab("Time (s)") +
  xlim(-50, 120) +
  theme_minimal() +
  theme(legend.position = "none")



# archive -----------------------------------------------------------------

lfout %>% subset(., predicted_subclass %in% c("L2/3_IT", "L5_IT", "L5_ET") & Cortical_area == "TCx") %>%
  mutate(x_bin = floor(t / 5) * 5) %>%
  group_by(x_bin, protocol, cell, predicted_subclass, Species) %>%
  summarise(y = mean(percent_change), .groups = "drop") %>%
  ggplot(., aes(x = x_bin, y = y, colour = cell)) +
  geom_line() +
  #facet_wrap( ~ predicted_subclass, nrow = 3, scales = "free") +
  facet_grid(cols = vars(Species), rows = vars(predicted_subclass)) +
  ylim(0, 300) +
  ylab("Time (s)") +
  xlim(-25, 120) +
  theme_minimal() +
  theme(legend.position = "none")





