#' Spike puff analysis to be nnested in the srt nnest
#'
#' @param x 5ht nnest
#'
#' @returns srt nnest
#' @export spikePuff5HT
#'
#' @examples
#' x = i
#'
#'
spikePuff5HT = function(x){

  #sh5HT = proj5HT::sheets5HT()
  MD = projHCT::sheets$MD
  #sMdta = sh5HT$sMdta

  sMdta = readODS::read_ods("den/serotonin/Data/selected_merged_metadata.ods", sheet = "Macaque")
  sMdta = sMdta[!is.na(sMdta$cell_name), ]
  sHdta = readODS::read_ods("den/serotonin/Data/selected_merged_metadata.ods", sheet = "Human")
  sHdta = sHdta[!is.na(sHdta$cell_name), ]

  mMD = bind_rows(sMdta,sHdta)



  if(is.character(x)){
  print(x)
  srt = loadHCT(x, tag = "-srt.rds")
  }else{
    x = srt$cell
  }

  iMD = MD[grep(x, MD$cell_name), ]
  iMerged = data.frame(mMD[grep(x, mMD$cell_name), ])
  iM0 = iMerged[which(grepl("spikeTTL", iMerged$stimulus_description)), ]

  if (nrow(iM0) > 1) {
    iM0 = iM0[grep(1, iM0$Rank), ]
  }

  if(!nrow(iM0)){
    message("No spikeTTL data for use")
    return(x)
  }

  if(is.na(iM0$keep)){
    message("check keep value. Keep is set NA. Need logical inferrence for analysis")
    return(x)
  }

  iM0 = iM0[iM0$keep, ]

  # If no protocols are going to be analyzed
  if (nrow(iM0) == 0) {
    return(x)
  }

  ttlTime = iM0$TTLtime

  if (is.na(ttlTime)) {
    ttlTime = 5
  }


  #pro = as.character(iM0[s, "stimulus_description"])
  pro = as.character(iM0$stimulus_description)
  spike_puff_sweeps = iM0$wrapped_sweeps
  spike_puff_sweeps = gsub("\\[|\\]", "", spike_puff_sweeps)
  spike_puff_sweeps = as.numeric(strsplit(spike_puff_sweeps, ",")[[1]])
  sN = spike_puff_sweeps + 1
  nwb = nphys::extractNWB(iMD$nwb_file, sweeps = sN)




  # For an unknown reason, sometimes the nwb returns the sweeps as not a
  if(is.null(names(nwb$sweeps))){
      nwb$sweeps = as.data.frame(nwb$sweeps)
  }

  sweep_names = names(nwb$sweeps)

  dur = sapply(names(nwb$sweeps), function(n) {
      sw = nwb$sweeps[[n]]
      #plot(sw,type = "l")
      sw = sw[which(!is.na(sw))]
      sr = nwb$sampling_rate[[n]]
      # if(sr == 25000){
      #   sr = 50000
      # }
      #plot(sw,type = "l")
      ret = tail(nphys::convertIndex(sw, sr), 1)
      # if(ret == 112){
      #   ret = 56
      # }else{
      #   ret = ret
      # }
    })

  sweep_peaks = sapply(names(nwb$sweeps), function(x) {
      #print(x)
      sweep = nwb$sweeps[[x]]
      #plot(sweep, type = 'l')
      # Drop NA values
      sweep = sweep[which(!is.na(sweep))]
      sr = nwb$sampling_rate[[x]]
      #print(sr)
      # if(sr == 25000){
      #   sr = 50000
      # }
      # run the find peaks function
      fp = nphys::APanalysis(sweep, rate = sr)

    })

  if(length(unlist(strsplit(pro, "\\+"))) != length(names(sweep_peaks))){
    names(sweep_peaks) = c(unlist(strsplit(pro, "\\+")), rep("trail", (length(names(sweep_peaks)) - length(unlist(strsplit(pro, "\\+"))))))
    }else{
    names(sweep_peaks) = unlist(strsplit(pro, "\\+"))
  }

  names(dur) = names(sweep_peaks)

  # Clear variables
    if (exists("dfBl")) {
      rm(dfBl)
    }
    if (exists("dfTTL")) {
      rm(dfTTL)
    }
    if (exists("dfTrail")) {
      rm(dfTrail)
    }
    maxpt = NULL
    maxpk = NULL

    if ("baseline" %in% names(sweep_peaks)) {
      fp = sweep_peaks[[grep("baseline", names(sweep_peaks))]]
      fp = as.data.frame(fp)
      #drop the last row... only necessary if using instRate instead of instRate_pre
      fp = fp[-nrow(fp), ]
      maxpk = as.numeric(dur[[grep("baseline", names(dur))]])
      epoch_count_time = fp$epoch_t - (maxpk + ttlTime)
      instRate = fp$instRate
      dfBl = data.frame(time = epoch_count_time, instRate = instRate)
      #plot(dfBl, ylim = c(0,25))
    }
    # Align the instantaneous rates of the test spikes so 0 lands on
    if ("spikeTTL" %in% names(sweep_peaks)) {
      fp = sweep_peaks[[grep("spikeTTL", names(sweep_peaks))]]
      fp = as.data.frame(fp)
      fp = fp[-nrow(fp), ]
      maxpt = as.numeric(dur[[grep("spikeTTL", names(dur))]])
      instRate = fp$instRate
      dfTTL = data.frame(time = fp$epoch_t - ttlTime, instRate = instRate)
      #plot(dfTTL)
    }
    # Align the instantaneous rates of the trail spikes
    if ("trail" %in% names(sweep_peaks)) {

      # If trail length is greater than 1 trail
      if(sum(names(sweep_peaks) == "trail") >= 2) {
        # Have to select the dataframes as a list instead of a dataframe
        fp = sweep_peaks[grep("trail", names(sweep_peaks))]
        fp = lapply(fp, function(l){
          l = l[-nrow(l), ]
        })
        fp = lapply(fp, as.data.frame)
        # Rep so that we can adjust according to how many trails we have
        lfp = unlist(sapply(1:length(names(fp)), function(l) {
          rep(l, nrow(fp[[l]]))
        }))
        #names(fp) = NULL
        fp = bind_rows(fp)

        # Use the trail count as a way to align the timings
        fp$trail_count = lfp
        # This adjust the trail x axis times for instances where a longer trail is present
        fp$epoch_t = sapply(1:nrow(fp), function(f) {
          fp$epoch_t[f] + as.numeric(dur[which(names(dur) == "trail")-1][fp$trail_count[f]]) *
            fp$trail_count[f] - ttlTime
        })
        fp$trail_count = NULL
      } else{
        fp = sweep_peaks[[grep("trail", names(sweep_peaks))]]
        fp = as.data.frame(fp)
        fp$epoch_t = fp$epoch_t + as.numeric(dur[which(names(dur) == "spikeTTL")]) -
          ttlTime
        fp = fp[-nrow(fp), ]
      }
      dfTrail = data.frame(time = fp$epoch_t,
                           instRate = fp$instRate)
    }

    tst = dfTTL

    if (exists("dfBl")) {
      tst = rbind.data.frame(dfBl, tst)
    }
    if (exists("dfTrail")) {
      tst = rbind.data.frame(tst, dfTrail)
    }

    # Remove any points that are greater than 5 standard deviations from the norm
    if(any(abs(tst$instRate - mean(tst$instRate, na.rm = TRUE)) > 5 * sd(tst$instRate, na.rm = TRUE))){
      tst = tst[-which(abs(tst$instRate - mean(tst$instRate, na.rm = TRUE)) > 5 * sd(tst$instRate, na.rm = TRUE)),]
    }

    tst$protocol = c(rep("baseline", length(which(tst$time < 0))), rep("stimulus", length(which(tst$time > 0))))

    blMean = mean(tst$instRate[which(tst$time < 0 & tst$time > -5)])
    tst$percent_change = (tst$instRate / blMean) * 100
    plot(tst$time, tst$percent_change, main = paste(x, iMerged$predicted[!is.na(iMerged$predicted)]), xlab = "time", ylab = "% Change")

    tst$cell_name = iMD$cell_name
    tst$predicted_subclass = iMD$predicted_subclass
    tst$assigned_subclass = iMerged$predicted[!is.na(iMerged$predicted)]
    tst$blockers = iM0$Blockers
    tst$ket = iM0$Ket
    tst$expCon = ifelse(is.na(iM0$Exp), "Standard_Puff", iM0$Exp)

    #tst$DFP = iMD$Depth_from_pia
    #tst$LIMs_depth = iMD$LIMS_depth
    #tst$assigned_depth = iMerged$depth[!is.na(iMerged$depth)]
    tst$response = iM0$Effect

    names(sN) = names(sweep_peaks)
    names(sweep_names) = names(sweep_peaks)



    out = list(
      cell = x,
      MD = MD,
      selected_MD = iMerged,
      pro_spikeTTL = iM0,
      sweep_names = sweep_names,
      spike_TTL_sweep_numbers = sN,
      ttlTime = ttlTime,
      sweep_duration = dur,
      sweep_peaks = sweep_peaks,
      spike_puff_output = tst,
      blMean = blMean
    )


    write.csv(tst, file.path(srt$rd, paste0(srt$cell, "-srtPuff.csv")))

    srt$dfs$spikeTTL = out
    saveRDS(srt, srt$files$nnest)

    print(x)

    return(srt)

}
