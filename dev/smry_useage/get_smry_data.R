smry = srt$exp$smry

smry$stimulus_description

library(nphys)
nwbPlot_sweep(srt$files$nwb, acquisition_name = smry$acquisition_names[44])

