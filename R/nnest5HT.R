#' Nnest the 5HT project data
#'
#' @param x
#' @param save_nnest
#'
#' @returns
#' @export
#'
#' @examples
#' x = i
#'
#'
nnest5HT = function(x, save_nnest = TRUE, update = TRUE){

  sheets = projHCT::sheets
  srt = list()
  carry = list()

  if (is.list(x)) {
    cell = x$cell
  } else if (is.character(x)) {
    cell = x
  }

  if(is.list(x) && update){
      if("dfs" %in% names(x)){
        carry$dfs = x$dfs
      }
      if("ana" %in% names(x)){
        carry$ana = x$ana
      }
  }

  md = sheets$MD[grep(cell, sheets$MD$cell_name),]
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
    cell = cell,
    wd = as.character(wd),
    rd = as.character(rd),
    md = md,
    project = "5HT",
    files = list(
      pxp = grep(cell, sheets$files$pxp_files, value = TRUE),
      nwb = grep(cell, sheets$files$nwb_files, value = TRUE),
      nnest = file.path(rd, paste0(cell, "-srt.rds")),
      hct_nnest = file.path(rd, paste0(cell, "-hct.rds")),
      md = file.path(rd, paste0(cell, projHCT::params$tags$md)),
      smry = file.path(rd, paste0(cell, projHCT::params$tags$smry))
    )
  )

  #### Summary dataframe
  if(file.exists(srt$files$smry)){
    srt$smry = read.csv(file.path(rd, paste0(cell, projHCT::params$tags$smry)))
  }

  if(update){
    srt = c(srt,carry)
  }

  if (file.exists(srt$files$md) && update) {
    write.csv(as.matrix(md),
              file = srt$files$md,
              row.names = FALSE)
  }

  if (!file.exists(srt$rd)) {
    message("creating new directory for cell analysis")
    dir.create(srt$rd)
    message("Writing Metadata")
    ## Write md.csv file point in nnest if it does not already exist
    if (!file.exists(srt$files$md)) {
      write.csv(as.matrix(md),
                file = srt$files$md,
                row.names = FALSE)
    }
  }


  if (save_nnest) {
    message(paste0("saving cell ", srt$cell,  "-srt.rds"))
    saveRDS(srt, file = srt$files$nnest)
  }
  return(srt)
}
