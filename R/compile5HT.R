#' Compile 5HT dataframes
#'
#' @param save_file
#'
#' @returns
#' @export
#'
#' @examples
compile5HT = function(save_file = TRUE){

  comp = list()

  lf = list.files("rookery", pattern = "-srtPuff.csv", recursive = TRUE, full.names = TRUE)

  lfout = lapply(lf, function(l) {
    print(l)
    #r = read.csv(l)
    r = readr::read_csv(l, show_col_types = FALSE)
    r$DFP = NULL
    r$LIMs_depth = NULL
    if(nrow(r)==0){
      return(NA)
    }
    return(r)
  }) %>% .[!is.na(.)] %>% bind_rows(.)

  lfout$...1 = NULL
  lfout$subclass_Corr = NULL

  if(save_file){
    write.csv(lfout, file = "data-raw/compile_srtPuff.csv", row.names = FALSE)
    }

  comp$srtPuff = lfout

  #return(lfout)

  lf = list.files("rookery", pattern = "-srtDrugWash.csv", recursive = TRUE, full.names = TRUE)

  lfout = lapply(lf, function(l) {
    print(l)
    #r = read.csv(l)
    r = readr::read_csv(l, show_col_types = FALSE)
    r$DFP = NULL
    r$LIMs_depth = NULL
    if(nrow(r)==0){
      return(NA)
    }
    return(r)
  }) %>% .[!is.na(.)] %>% bind_rows(.)

  lfout$...1 = NULL
  lfout$subclass_Corr = NULL

  if(save_file){
    write.csv(lfout, file = "data-raw/compile_srtDrugWash.csv", row.names = FALSE)
  }

  comp$srtDrugWash = lfout

  return(comp)


}
