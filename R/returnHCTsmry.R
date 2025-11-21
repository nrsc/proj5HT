#' Return HCT smry
#'
#' @param x
#'
#' @returns
#' @export
#'
#' @examples
returnHCTsmry = function(x){
  hct = loadHCT(x)

  smry = hct$exp$smry

  return(smry)

}
