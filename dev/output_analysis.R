### Script to present output from Brians python analysis
MD = projHCT::sheets$MD
cells = list.files("den/serotonin/output/")
cell_json = list.files("den/serotonin/output/", recursive = TRUE, full.names = TRUE, pattern = ".json")
merged_metadata = read.csv("den/serotonin/Data/merged_metadata.csv")
output_cellsMD = MD[MD$cell_name %in% cells,]
length(output_cellsMD$cell_name)

#which
unique(merged_metadata$cell_name) %in% cells
output_cellsMD$cell_name %in% cells
length(output_cellsMD$cell_name %in% cells)

#### Reading the json ####
# my_list = jsonlite::fromJSON(cell_json[1])
# my_list = jsonlite::fromJSON(cell_json[3])
# plot(x = my_list$spike_puff0$spike_times, y = my_list$spike_puff0$Instantaneous_frequency)

table(output_cellsMD$predicted_subclass)
cells5ET = output_cellsMD[which(output_cellsMD$predicted_subclass == "L5_ET"),]

srt$exp$smry$stimulus_description
i = cells5ET$cell_name[5]



for(i in cells5ET$cell_name[3:length(cells5ET$cell_name)]){

  iMD = MD[grep(i, MD$cell_name),]
  srt = loadHCT(i)
  smry = srt$exp$smry
  iMerged = merged_metadata[grep(i, merged_metadata$cell_name),]
  jfile = cell_json[grep(i, cell_json)]
  jlist = jsonlite::fromJSON(jfile)
  sweeps_analyzed = jlist[[1]]$sweeps_analyzed

  if("baseline_spike_puff" %in% iMerged$stimulus_description){


  baseline_spike_puff = iMerged[which(iMerged$stimulus_description == "baseline_spike_puff"), "wrapped_sweeps"]
  baseline_spike_puff = gsub("\\[|\\]", "", baseline_spike_puff)
  baseline_spike_puff = as.numeric(strsplit(baseline_spike_puff, ",")[[1]])
  n = baseline_spike_puff[1]

  spike_puff_sweeps = iMerged[which(iMerged$stimulus_description == "spike_puff"), "wrapped_sweeps"]
  spike_puff_sweeps = gsub("\\[|\\]", "", spike_puff_sweeps)
  spike_puff_sweeps = as.numeric(strsplit(spike_puff_sweeps, ",")[[1]])
  n2 = spike_puff_sweeps[1]


  sweep = extractHCT(srt$files$nwb, sweeps = n+1, slim = TRUE)

  smry[n+1,]

    plot(sweep, type = "l")

    fp = nphys::APanalysis(sweep, rate = srt$exp$smry$sampling_rate[n+3])
    fp = projHCT::APanalysisHCT(sweep, rate = srt$exp$smry$sampling_rate[n])
    fp = as.data.frame(fp)
    fp$epoch_t
    bin_width = 5


    out1 = fp %>% mutate(bin = cut(epoch_t, breaks = seq(0, max(fp$epoch_t), by = bin_width), right = FALSE)) %>%
      group_by(bin) %>%
      summarise(count = n(),
                rate = count / bin_width)
    out1 = out1[(nrow(out1)-5):nrow(out1),]

    tstSweep = extractHCT(srt$files$nwb, sweeps = n2+1, slim = TRUE)
    plot(tstSweep, type = "l")
    fp2 = nphys::APanalysis(tstSweep, rate = srt$exp$smry$sampling_rate[n2+1])
    fp2 = as.data.frame(fp2)
    #bins2 <- cut(fp2$epoch_t, breaks = breaks, right = FALSE)
    #table(bins2)
    out2 = fp2 %>% mutate(bin = cut(
      epoch_t,
      breaks = seq(0, 55, by = bin_width),
      right = FALSE
    )) %>%
      group_by(bin) %>%
      summarise(count = n(), rate = count / bin_width)

    # Drop last row
    out2 = out2[-nrow(out2), ]



    #}
    # Stimulus sweep
    #plot(tstSweep, type = "l")
    # Bin averages of baseline sweep
    pCh = as.numeric(out1[nrow(out1), "rate"])
    df = rbind(out1, out2)

    df$pro = c(rep("baseline", nrow(out1)), rep("stimulus", nrow(out2)))
    df$percent_change = df$rate/pCh*100
    df$bin = NULL
    df$count = NULL
    df = as.data.frame(df)
    df$time = seq(-30,50, by = 5)
    gg = ggplot(df, aes(x = time, y = percent_change, colour = pro)) +
      geom_point() +
      ylim(0, max(df$percent_change)) +
      theme_minimal()

    plot(gg)
    #plotVoid(tst)
    #plotVoid(tst2)

    df$cell_name = srt$cell
    df$predicted_subclass = srt$md$predicted_subclass
}

    readline("write csv")
    write.csv(df, file.path(srt$rd, paste0(srt$cell, "-srt50sPuff.csv")))

  }

lf = list.files("rookery", pattern = "-srt50sPuff.csv", recursive = TRUE, full.names = TRUE)


df = read.csv(lf[3])

df

gg = ggplot(df, aes(x = time, y = percent_change, colour = pro)) +
  geom_point() +
  ylim(0, max(df$percent_change)) +
  ggtitle(unique(df$cell_name)) +
  theme_minimal()

plot(gg)



#### Build artifacts
# Baseline Sweep
#n = iMerged[which(iMerged$stimulus_description == "baseline_spike_puff"), "wrapped_sweeps"]
#n = as.numeric(gsub("\\[|\\]", "", n))
# if(length(n) > 1){
#   print(jlist[[1]]$sweeps_analyzed)
#   if(n[1]+1 == swanal[1]){
#     n = n[1]
#   }
# }
# Stimulus sweeps
#iMerged[which(iMerged$stimulus_description == "spike_puff"), "wrapped_sweeps"]

# if(jlist[[1]]$sweeps_analyzed[1] == n+1){
#   n2 = jlist[[1]]$sweeps_analyzed[1]
# }else{
#   print(n)
#   print(jlist[[1]]$sweeps_analyzed)
#   n2 = readline("select which puff sweep to analyze")
#if(jlist[[1]]$sweeps_analyzed-1 == n){
#n2 = iMerged[which(iMerged$stimulus_description == "spike_puff"), "wrapped_sweeps"]
#n2 = jlist[[1]]$sweeps_analyzed
#n2 = as.numeric(gsub("\\[|\\]", "", n2))


plot(x = jlist$spike_puff0$spike_times, y = jlist$spike_puff0$Instantaneous_frequency)
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


if(smry$stimulus_description[preTTL] == smry$stimulus_description[rSweeps[1]])



### For future analysis
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


