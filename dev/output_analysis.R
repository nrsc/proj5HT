### Script to present output from Brians python analysis
MD = projHCT::sheets$MD
cells = list.files("den/serotonin/output/")
cell_json = list.files("den/serotonin/output/", recursive = TRUE, full.names = TRUE, pattern = ".json")
merged_metadata = read.csv("den/serotonin/Data/merged_metadata.csv")
output_cellsMD = MD[MD$cell_name %in% cells,]

length(output_cellsMD$cell_name)

unique(merged_metadata$cell_name) %in% cells

output_cellsMD$cell_name %in% cells
length(output_cellsMD$cell_name %in% cells)

### Not sure why this cell is not showing up in my projHCT::MD. The JNB has it labeled H22.26.429.12.57.03, but it is in the MD file as H22.26.429.12.57.01.03
#grep(cells[25], MD$cell_name) # Fixed by running sheetsHCT()

msmry = projHCT::sheets$dfs$massSmry



#### Cells with 50s puff protocols ####
tst = msmry[grep("50spuff", msmry$stimulus_description),]
unique(tst$cell_name)
length(unique(tst$cell_name))

#### mass summary data for cells already analyzed
tst2 = msmry[msmry$cell_name %in% cells,]
unique(tst2$cell_name)
unique(tst2$stimulus_description)


#### MD descriptions for all 50s puff cells ####
MDall = MD[MD$cell_name %in% unique(tst$cell_name),]
table(MDall$Species)

table(gsub(" ", "_", MDall$subclass_Corr))
table(MDall$predicted_subclass)
length(which(is.na(MDall$predicted_subclass)))

# Cells not already analyzed
for_analysis = unique(tst$cell_name)[!unique(tst$cell_name) %in% unique(tst2$cell_name)]
MDanalysis = MD[MD$cell_name %in% for_analysis,]




#### Reading the json ####
my_list = jsonlite::fromJSON(cell_json[1])
my_list = jsonlite::fromJSON(cell_json[3])
plot(x = my_list$spike_puff0$spike_times, y = my_list$spike_puff0$Instantaneous_frequency)
i = cells[3]


for(i in cells){

  iMD = MD[grep(i, MD$cell_name),]
  srt = loadHCT(i)
  iMerged = merged_metadata[grep(i, merged_metadata$cell_name),]
  jfile = cell_json[grep(i, cell_json)]
  jlist = jsonlite::fromJSON(jfile)

  if("baseline_spike_puff" %in% iMerged$stimulus_description){

    # Baseline Sweep
    n = iMerged[which(iMerged$stimulus_description == "baseline_spike_puff"), "wrapped_sweeps"]
    n = as.numeric(gsub("\\[|\\]", "", n))
    tst = extractHCT(srt$files$nwb, sweeps = n+1, slim = TRUE)
    #plot(tst, type = "l")
    fp = nphys::APanalysis(tst, rate = srt$exp$smry$sampling_rate[n+1])
    fp = as.data.frame(fp)

    # Bin averages of baseline sweep
    out1 = fp %>% mutate(bin = cut(epoch_t, breaks = seq(0, max(epoch_t), by = bin_width), right = FALSE)) %>%
      group_by(bin) %>%
      summarise(count = n(),
                rate = count / bin_width)
    out1 = out1[-nrow(out1),]

    # Stimulus sweep
    n2 = iMerged[which(iMerged$stimulus_description == "spike_puff"), "wrapped_sweeps"]
    n2 = as.numeric(gsub("\\[|\\]", "", n2))
    tst2 = extractHCT(srt$files$nwb, sweeps = n2+1, slim = TRUE)
    #plot(tst2, type = "l")

    fp2 = nphys::APanalysis(tst2, rate = srt$exp$smry$sampling_rate[n2+1])
    fp2 = as.data.frame(fp2)
    bins2 <- cut(fp2$epoch_t, breaks = breaks, right = FALSE)

    table(bins2)

    out2 = fp2 %>% mutate(bin = cut(epoch_t, breaks = seq(0, 55, by = bin_width), right = FALSE)) %>%
      group_by(bin) %>%
      summarise(count = n(),
                rate = count / bin_width)
    out2 = out2[-nrow(out2),]


    pCh = as.numeric(out1[nrow(out1), "rate"])
    df = rbind(out1[6:11,], out2)
    df$pro = c(rep("baseline", 7), rep("stimulus", 10))

    df$percent_change = df$rate/pCh*100

    df$bin = NULL
    df$count = NULL

    df = as.data.frame(df)
    df$time = seq(-30,50, by = 5)

    ggplot(df, aes(x = time, y = percent_change, colour = pro)) +
      geom_point() +
      ylim(0, max(df$percent_change)) +
      theme_minimal()

    #plotVoid(tst)
    #plotVoid(tst2)

    comparison <- data.frame(
      bin = 1:n_bins,
      stim_rate = stim_rate,
      baseline_rate = rep(baseline_rate, n_bins)
    )
  }

  df = data.frame(time = fp$epoch_t)
  bin_width = 5

  df %>% mutate(bin = cut(time, breaks = seq(0, max(time)+bin_width, by = bin_width), right = FALSE)) %>%
    group_by(bin) %>%
    summarise(count = n(),
              rate = count / bin_width)

  smry = srt$exp$smry
  smry[jlist$spike_puff0$sweeps_analyzed+1,]

  rSweeps = jlist$spike_puff0$sweeps_analyzed+1
  preTTL = rSweeps[1]-1
  grep(smry[my_list$spike_puff0$sweeps_analyzed+1, "stimulus_description"][1], smry$stimulus_description)


  if(smry$stimulus_description[preTTL] == smry$stimulus_description[rSweeps[1]]){

  }

  plot(x = jlist$spike_puff0$spike_times, y = jlist$spike_puff0$Instantaneous_frequency)

}


