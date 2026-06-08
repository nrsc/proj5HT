NaAvail5HT = function(x, show_plot = TRUE, return_dfs = TRUE){


np <- reticulate::import("numpy")
tool <- reticulate::import_from_path("sodium_tool", "py")
signal <- reticulate::import("scipy.signal")


if(is.character(x)){
  print(x)
  srt = load5HT(x, tag = "-srt.rds", rehydrate = FALSE)
}else{
  srt = x
}


blSweep_name = srt$dfs$spikeTTL$sweep_names[[grep("baseline", names(srt$dfs$spikeTTL$sweep_name))]]
blpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("baseline", names(srt$dfs$spikeTTL$sweep_peaks))]]
blsweep = nphys::extractNWB(srt$md$nwb_file, sweep = blSweep_name, slim = TRUE)
blrate = srt$smry$sampling_rate[[grep(blSweep_name, srt$smry$acquisition_names)]]


sweep_name = srt$dfs$spikeTTL$sweep_names[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_name))]]
TTLpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_peaks))]]
sweep = nphys::extractNWB(srt$files$nwb, sweep = sweep_name, slim = TRUE)
rate = srt$smry$sampling_rate[[grep(sweep_name, srt$smry$acquisition_names)]]

trailSweep_names = srt$dfs$spikeTTL$sweep_names[[grep("trail", names(srt$dfs$spikeTTL$sweep_name))]][1]
trailpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("trail", names(srt$dfs$spikeTTL$sweep_peaks))]]
trailsweep = nphys::extractNWB(srt$md$nwb_file, sweep = trailSweep_names, slim = TRUE)
trailrate = srt$smry$sampling_rate[[grep(trailSweep_names, srt$smry$acquisition_names)]]

#str(trailsweep)

# Make sure both sweeps have the same sampling rate (resample if needed)
#


if (rate != trailrate) {
  message("Resampling trailing sweep to match main sweep rate...")

  trailsweep = trailsweep[!is.na(trailsweep)]
  trailsweep_vec <- as.numeric(trailsweep)
  n_old <- length(trailsweep_vec)
  scale_factor <- rate / trailrate
  n_new <- round(as.numeric(n_old) * scale_factor)  # keep as double, not integer
  if (!is.finite(n_new) || n_new <= 0) {
    stop(sprintf("Invalid n_new: n_old=%g, rate=%g, trailrate=%g", n_old, rate, trailrate))
  }

  # Pass the Python integer safely
  trailsweep <- signal$resample(trailsweep_vec, as.integer(n_new))
  message(sprintf("Resampled trailing sweep: %.0f → %d samples", n_old, n_new))
}
cat("Rate:", rate, " Trailrate:", trailrate, "\n")
cat("nrow(trailsweep):", NROW(trailsweep), "\n")

drug_combined <- np$concatenate(list(np$array(as.numeric(sweep)), np$array(as.numeric(trailsweep))))

length(blsweep) > 0
length(sweep) > 0
length(trailsweep) > 0
blrate > 0
rate > 0
trailrate > 0
length(sweep) + length(trailsweep) == length(drug_combined)

# Concatenate them seamlessly

# #### Fit just drug ####

#fit_drug <- fit_params_from_dr(

# fit_full <- fit_params_from_baseline(
#     base = np$ravel(np$array(blsweep)),
#   #v_drug = np$ravel(np$array(sweep)),
#   fs = blrate,
#   #fs_drug = rate,
#   spike_thresh = -20L
# )

# res_pair <- infer_availability_continuous_pair(
#   v_base = np$ravel(np$array(blsweep)),
#   v_drug = np$ravel(np$array(sweep)),
#   fitted = fit_full,
#   fs_base = blrate,
#   fs_drug = rate,
#   spike_thresh = -20L
# )


fit_drug <- tool$fit_params_from_drug_segment(
  v_drug = np$array(sweep),
  fs_drug = rate,
  spike_thresh = -20L
)

#drug_plus_trail <- np$concatenate(list(np$array(sweep), np$array(trailsweep)))

res_pair <- tool$infer_availability_continuous_pair(
  v_base = np$array(blsweep),
  v_drug = np$array(drug_combined),
  #v_drug = np$array(sweep),
  fitted = fit_drug,          # from fit_params_from_baseline or FI/VI
  fs_base = blrate,
  fs_drug = rate,
  spike_thresh = -20L
)


df_b <- data.frame(time = res_pair$time_base, a = res_pair$a_cont_base, cond="Baseline")
df_d <- data.frame(time = res_pair$time_drug, a = res_pair$a_cont_drug, cond="Drug")


#tail(df_d$time,100)

#max(df_d$time)

# assume spike_puff_output$time has 0 = drug application
t_shift <- min(srt$dfs$spikeTTL$spike_puff_output$time, na.rm = TRUE)

#--- 2. Shift sodium availability times so baseline aligns negative ----
df_b$time <- df_b$time + t_shift
#df_b$time <- srt$dfs$NaAv$df_b$time + t_shift
df_d$time <- df_d$time + t_shift + max(df_b$time) - min(df_b$time, na.rm = TRUE)
#df_d$time <- srt$dfs$NaAv$df_d$time + t_shift + max(srt$dfs$NaAv$df_b$time) - min(srt$dfs$NaAv$df_b$time, na.rm = TRUE)

#--- 3. Combine sodium availability ----
df_na <- rbind(df_b, df_d)



#range(df_b$a)


#range(srt$dfs$spikeTTL$spike_puff_output$time)

#--- 4. Prepare %-change points, rescaled for 0–1 overlay ----
df_spike <- srt$dfs$spikeTTL$spike_puff_output
names(df_spike)[names(df_spike) == "percent_change"] <- "percent_change"
df_spike$percent_scaled <- scales::rescale(df_spike$percent_change, to = c(0, 1))

if(show_plot){

#--- 5. Plot ----
gg = ggplot() +
  # sodium availability background
  geom_line(
    data = df_na,
    aes(x = time, y = a, color = cond),
    linewidth = 0.6, alpha = 0.7
  ) +
  # spike %-change points (rescaled)
  geom_point(
    data = df_spike,
    aes(x = time, y = percent_scaled),
    shape = 21, fill = "black", color = "white",
    size = 2, alpha = 0.9
  ) +
  # vertical marker for drug onset
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  annotate("text", x = 0, y = 1.02, label = "Drug onset", hjust = -0.1, size = 4) +
  labs(
    title = paste(
      "Continuous Sodium Availability and Spike % Change —",
      srt$dfs$spikeTTL$selected_MD$cell_name
    ),
    x = "Time (s, 0 = drug onset)",
    y = "Sodium Availability a(t)"
  ) +
  scale_y_continuous(
    name = "Sodium Availability a(t)",
    sec.axis = sec_axis(
      trans = ~ . * (max(df_spike$percent_change, na.rm = TRUE) -
                       min(df_spike$percent_change, na.rm = TRUE)) +
        min(df_spike$percent_change, na.rm = TRUE),
      name = "% Change"
    )
  ) +
  scale_color_manual(values = c("Baseline" = "salmon", "Drug" = "cyan3")) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

plot(gg)

}

if(return_dfs){

dfs = list(
  fit_drug = fit_drug,
  df_b = df_b,
  df_d = df_d,
  df_na = df_na,
  df_spike = df_spike
)

return(dfs)
}

return(NULL)

}


