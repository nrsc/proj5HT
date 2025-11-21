library(reticulate)
library(ggplot2)

np <- import("numpy")
tool <- import_from_path("sodium_tool", "py")

srt = loadHCT("Q20.26.026.2A.05.02", tag = "-srt.rds")

blSweep_name = srt$dfs$spikeTTL$sweep_names[[grep("baseline", names(srt$dfs$spikeTTL$sweep_name))]]
blpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("baseline", names(srt$dfs$spikeTTL$sweep_peaks))]]
blsweep = nphys::extractNWB(srt$md$nwb_file, sweep = blSweep_name, slim = TRUE)
blrate = srt$smry$sampling_rate[[grep(blSweep_name, srt$smry$acquisition_names)]]

sweep_name = srt$dfs$spikeTTL$sweep_names[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_name))]]
TTLpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_peaks))]]
sweep = nphys::extractNWB(srt$files$nwb, sweep = sweep_name, slim = TRUE)
rate = srt$smry$sampling_rate[[grep(sweep_name, srt$smry$acquisition_names)]]

fit_base <- tool$fit_params_from_baseline(np$array(blsweep), fs = fs, spike_thresh = -20L)

res_pair <- tool$infer_availability_continuous_pair(
  v_base = np$array(blsweep),
  v_drug = np$array(sweep),
  fitted = fit_base,          # from fit_params_from_baseline or FI/VI
  fs_base = blrate,
  fs_drug = rate,
  spike_thresh = -20L
)

# Continuous a(t)
df_b <- data.frame(time = res_pair$time_base, a = res_pair$a_cont_base, cond="Baseline")
df_d <- data.frame(time = res_pair$time_drug, a = res_pair$a_cont_drug, cond="Drug")

ggplot(rbind(df_b, df_d), aes(time, a, color = cond)) +
  geom_line() +
  labs(title="Continuous Sodium Availability", x="Time (s)", y="a(t)") +
  theme_minimal()

# Fit from baseline

# Plot

# df_b <- data.frame(time = res_pair$time_base, a = res_pair$a_cont_base, cond="Baseline")
# df_d <- data.frame(time = res_pair$time_drug, a = res_pair$a_cont_drug, cond="Drug")
# ggplot(rbind(df_b,df_d), aes(time,a,color=cond)) +
#   geom_line() + theme_minimal() +
#   labs(title="Continuous Sodium Availability", y="a(t)", x="Time (s)")



# #### Fit just drug ####
# fit_drug <- tool$fit_params_from_drug_segment(
#   v_drug = np$array(sweep),
#   fs_drug = rate,
#   spike_thresh = -20L
# )
#
# res_pair <- tool$infer_availability_continuous_pair(
#   v_base = np$array(blsweep),
#   v_drug = np$array(sweep),
#   fitted = fit_drug,          # from fit_params_from_baseline or FI/VI
#   fs_base = blrate,
#   fs_drug = rate,
#   spike_thresh = -20L
# )
#
# df_b <- data.frame(time = res_pair$time_base, a = res_pair$a_cont_base, cond="Baseline")
# df_d <- data.frame(time = res_pair$time_drug, a = res_pair$a_cont_drug, cond="Drug")
# ggplot(rbind(df_b,df_d), aes(time,a,color=cond)) +
#   geom_line() + theme_minimal() +
#   labs(title="Continuous Sodium Availability", y="a(t)", x="Time (s)")


