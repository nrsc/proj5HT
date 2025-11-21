analysis5HT_dvdt5HT = function(x, t0 = 25){

  srt = loadHCT(x, tag = "-srt.rds")

  blSweep_name = srt$dfs$spikeTTL$sweep_names[[grep("baseline", names(srt$dfs$spikeTTL$sweep_name))]]
  blpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("baseline", names(srt$dfs$spikeTTL$sweep_peaks))]]
  blsweep = extractHCT(srt$md$nwb_file, sweep = blSweep_name, slim = TRUE)


  sweep_name = srt$dfs$spikeTTL$sweep_names[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_name))]]
  TTLpeaks = srt$dfs$spikeTTL$sweep_peaks[[grep("spikeTTL", names(srt$dfs$spikeTTL$sweep_peaks))]]
  sweep = extractHCT(srt$files$nwb, sweep = sweep_name, slim = TRUE)

  rate = srt$smry$sampling_rate[[grep(sweep_name, srt$smry$acquisition_names)]]
  fp = as.data.frame(TTLpeaks)
  fp = fp[-nrow(fp),]

  preAPs = sweep[head(fp$apONi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime - 5],1):tail(fp$apOFFi[fp$epoch_t < srt$dfs$spikeTTL$ttlTime],1)]
  postAPs = sweep[head(fp$apONi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime],1):head(fp$apOFFi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime + 5],1)]

  #t0 = 20#fp$epoch_t[which(fp$instRate == max(fp$instRate))]
  tLess = t0-2.5
  tPlus = t0+2.5

  maxAPs = sweep[head(fp$apONi[fp$epoch_t > tLess],1):head(fp$apOFFi[fp$epoch_t > tPlus],1)]

  par(mfrow = c(1,3))

  trace1 = dvdt_mV(preAPs, rate = rate, return = TRUE, plot = FALSE)

  trace2 = dvdt_mV(postAPs, rate = rate, return = TRUE, plot = FALSE)

  traceMax = dvdt_mV(maxAPs, rate = rate, return = TRUE, plot = FALSE)

  df1 = data.frame(trace1$dsignal)
  df2 = data.frame(trace2$dsignal)
  df3 = data.frame(traceMax$dsignal)

  df = bind_rows(df1,df2,df3)

  df0$ref = c(rep("baseline", nrow(df1)),rep("TTL", nrow(df2)),rep("max", nrow(df3)))

  gg0 = ggplot(df0, aes(x = L, y = V2, colour = ref)) +
    geom_path()

  par(mfrow = c(1,3))

  firstBLAP = sweep[head(fp$apONi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime - 5],1):head(fp$apOFFi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime - 5],1)]

  TTLAP = sweep[head(fp$apONi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime],1):head(fp$apOFFi[fp$epoch_t > srt$dfs$spikeTTL$ttlTime],1)]

  maxrateAP = sweep[head(fp$apONi[fp$epoch_t >= t0],1):head(fp$apOFFi[fp$epoch_t >= t0],1)]

  trace1 = dvdt_mV(firstBLAP, rate = rate, return = TRUE)

  trace2 = dvdt_mV(TTLAP, rate = rate, return = TRUE)

  traceMax = dvdt_mV(maxrateAP, rate = rate, return = TRUE)

  df1 = data.frame(trace1$dsignal)
  df2 = data.frame(trace2$dsignal)
  df3 = data.frame(traceMax$dsignal)

  df0 = bind_rows(data.frame(trace1$dsignal),data.frame(trace2$dsignal),data.frame(traceMax$dsignal))

  df0$ref = c(rep("baseline", nrow(df1)),rep("TTL", nrow(df2)),rep("max", nrow(df3)))



  ggSingleAP = ggplot(df0, aes(x = L, y = V2, colour = ref)) +
    geom_path()





  # If we remember from the TTL section of the sweep there is a delay between when the TTL goes on and the 5HT hits

#plot(fp$dvdtMax)

#plot(fp$dvdtMin)




}
