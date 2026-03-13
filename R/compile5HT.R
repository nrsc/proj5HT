#' Compile 5HT dataframes
#'
#' @param save_file
#'
#' @returns
#' @export
#'
#' @examples
compile5HT = function(save_file = TRUE){

  library(stringr)

  MD = projHCT::sheets$MD

  comp = list()

  lf = list.files("rookery", pattern = "-srtPuff.csv", recursive = TRUE, full.names = TRUE)

  lfout = lapply(lf, function(l) {
    print(l)
    #r = read.csv(l)
    r = readr::read_csv(l, show_col_types = FALSE)
    #r$DFP = NULL
    #r$LIMs_depth = NULL
    if(nrow(r)==0){
      return(NA)
    }
    return(r)
  }) %>% .[!is.na(.)] %>% bind_rows(.)

  #unique(lfout$DFP)

  ### Dataframe modifications
  lfout$...1 = NULL
  #lfout$subclass_Corr = NULL
  lfout$bath = gsub("Way635", "WAY635", lfout$bath)
  #unique(lfout$assigned_subclass)

  outMD = MD[MD$cell_name %in% unique(lfout$cell_name),]


  #outMD$LIMS_depth = as.numeric(outMD$LIMS_depth)
  #MD$assigned_depth = #if(is.na(outMD$LIMS_depth))

  comp$srtPuff = lfout

  comp$srtPuff <- comp$srtPuff %>%
    left_join(
      MD %>% select(cell_name, Date),
      by = "cell_name"
    )

  if(save_file){
    write.csv(comp$srtPuff, file = "data-raw/compile_srtPuff.csv", row.names = FALSE)
    }


  lf_dw = list.files("rookery", pattern = "-srtDrugWash.csv", recursive = TRUE, full.names = TRUE)

  lfout_dw = lapply(lf_dw, function(l) {
    print(l)
    #r = read.csv(l)
    r = readr::read_csv(l, show_col_types = FALSE)
    #r$DFP = NULL
    #r$LIMs_depth = NULL
    if(nrow(r)==0){
      return(NA)
    }
    return(r)
  }) %>% .[!is.na(.)] %>% bind_rows(.)


  ### Dataframe modifications
  lfout_dw$...1 = NULL
  #lfout$subclass_Corr = NULL
  lfout_dw$bath = gsub("Way635", "WAY635", lfout_dw$bath)


  comp$srtDrugWash = lfout_dw
  comp$srtDrugWash <- comp$srtDrugWash %>%
    left_join(
      MD %>% select(cell_name, Date),
      by = "cell_name"
    )
  if(save_file){
    write.csv(comp$srtDrugWash, file = "data-raw/compile_srtDrugWash.csv", row.names = FALSE)
  }

  saveRDS(comp, "data-raw/compiled5HT.rds")

  return(comp)


}
