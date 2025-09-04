nnest5HT = function(x){

  srt = list()

  if (is.list(x)) {
    cell = x$cell
  } else if (is.character(x)) {
    cell = x
  }

  #cell = "H22.03.304.11.03.01.01"

  sheets = projHCT::sheets
  #MD = sheets$MD
  #cMD = sheets$cleanMD

  md = sheets$cleanMD[grep(cell, sheets$cleanMD$cell_name),]
  if(nrow(md)==0){
    md = sheets$MD[grep(cell, sheets$MD$cell_name),]

    md = as.data.frame(md)
    if(md$User == "Nikolai"){
      srt = niks_nnestHCT(cell)
      return(srt)
    }
    if(file.exists(md$pxp_file)){
      wd = nphys::wrkD(md$pxp_file)
    }else if(file.exists(md$nwb_file)){
      wd = nphys::wrkD(md$nwb_file)
    }else{
      return(NULL)
    }
  }else{

    #md = sheets$cleanMD[grep(cell, sheets$MD$cell_name),]
    md = as.data.frame(md)
    wd = nphys::wrkD(md$pxp_file)
  }

  rd = file.path(projHCT::params$rookery, md$cell_name)

  srt = list(
    project = "5HT",
    cell = cell,
    wd = as.character(wd),
    rd = as.character(rd),
    md = md,
    files = list(
      pxp = grep(cell, sheets$files$pxp_files, value = TRUE),
      nwb = grep(cell, sheets$files$nwb_files, value = TRUE),
      nnest = file.path(rd, paste0(cell, "-srt.rds")),
      md = file.path(rd, paste0(cell, projHCT::params$tags$md)),
      smry = file.path(rd, paste0(cell, projHCT::params$tags$smry))
    )
  )

  #### Experimental parameters here
  #### Summary dataframe here



  if (!file.exists(srt$rd)) {
    message("creating new directory for cell analysis")
    dir.create(srt$rd)
    message("Writing Metadata")
    ## Write mdHCT.csv file point in nnest if it does not already exist
    if (!file.exists(srt$files$md)) {
      write.csv(as.matrix(md),
                file = srt$files$md,
                row.names = FALSE)
    }
  }

  if (!file.exists(srt$files$md)) {
    write.csv(as.matrix(md),
              file = srt$files$md,
              row.names = FALSE)
  }


  return(srt)




}
