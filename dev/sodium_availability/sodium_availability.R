srt = loadHCT("Q20.26.026.2A.05.02", tag = "-srt.rds")

blSweep_name = srt$dfs$spikeTTL$sweep_names[[grep("baseline", names(srt$dfs$spikeTTL$sweep_name))]]
blpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("baseline", names(srt$dfs$spikeTTL$sweep_peaks))]]
blsweep = nphys::extractNWB(srt$md$nwb_file, sweep = blSweep_name, slim = TRUE)

sweep_name = srt$dfs$spikeTTL$sweep_names[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_name))]]
TTLpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_peaks))]]
sweep = nphys::extractNWB(srt$files$nwb, sweep = sweep_name, slim = TRUE)

rate = srt$smry$sampling_rate[[grep(sweep_name, srt$smry$acquisition_names)]]
fp = as.data.frame(TTLpeaks)



#### Run analysis

np <- import("numpy")
voltage <- np$array(as.numeric(sweep))
#voltage = as.numeric(sweep)
fs = rate

reticulate::source_python('~/proj5HT/py/sodium_availability.py')

# detect spikes
spike_idx <- detect_spikes(voltage, fs, thresh=-20)

# extract spike features
features_list <- extract_spike_features(voltage, spike_idx, fs)

# infer sodium availability
w_vec <- c(0.01, -0.5, 0.001, -0.1)
w0 <- -2
tau_rec <- 0.05

a_trace <- infer_sodium_availability(features_list, w_vec, w0, tau_rec)


plot(sapply(features_list, function(f) f$ISI), a_trace, type="b",
     xlab="ISI", ylab="Inferred sodium availability")

plot(1:length(a_trace), a_trace, type="b",
     xlab="Spike number", ylab="Inferred sodium availability")
