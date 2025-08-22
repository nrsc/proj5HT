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


i = cells[1]

for(i in cells){

  iMD = MD[grep(i, MD$cell_name),]
  srt = loadHCT(i)

  srt$f$nwb

  jfile = cell_json[grep(i, cell_json)]
  jlist = jsonlite::fromJSON(jfile)

  smry = srt$exp$smry
  smry[jlist$spike_puff0$sweeps_analyzed+1,]

  rSweeps = jlist$spike_puff0$sweeps_analyzed+1
  preTTL = rSweeps[1]-1
  grep(smry[my_list$spike_puff0$sweeps_analyzed+1, "stimulus_description"][1], smry$stimulus_description)


  if(smry$stimulus_description[preTTL] == smry$stimulus_description[rSweeps[1]]){

  }

  plot(x = jlist$spike_puff0$spike_times, y = jlist$spike_puff0$Instantaneous_frequency)

}

msry = projHCT::sheets$dfs$massSmry

unique(msry$stimulus_description)

puffPros = msmry[grep("50spuff", msmry$stimulus_description),]

length(unique(puffPros$cell_name))



