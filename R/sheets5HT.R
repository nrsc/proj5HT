#' Sheets critical for 5HT project analysis functions
#'
#'
#' @returns
#' @export
#'
#' @examples
sheets5HT = function(){

  sMdta = readODS::read_ods("~/proj5HT/den/serotonin/Data/selected_merged_metadata.ods", sheet = "Macaque")
  sMdta = sMdta[!is.na(sMdta$cell_name), ]






  sHdta = readODS::read_ods("~/proj5HT/den/serotonin/Data/selected_merged_metadata.ods", sheet = "Human")
  sHdta = sHdta[!is.na(sHdta$cell_name), ]





  sh5HT = list(
    mUI = sMdta,
    hUI = sHdta
    )

  nnested_cells = list.files(
    "rookery",
    pattern = "-srt.rds",
    recursive = TRUE,
    full.names = TRUE
  ) %>% nphys::rmFlag(.) %>% sapply(., nphys::fileD) %>% as.character(.)

  sh5HT$cells = cells
  sh5HT$qCells = cells[grep("Q", cells)]
  sh5HT$hCells = cells[grep("h", cells)]

  if(grepl("proj5HT", rprojroot::find_rstudio_root_file())) {
    if (!identical(proj5HT::sh5HT, sh5HT)) {

      message("There are internal project data does not match updated sheets data")
      if (askYesNo("Would you like to update the package data (this is recommended")) {
        # Update sheets data to project
        message("updating sheets5HT built in data")

        usethis::use_data(sh5HT, overwrite = TRUE)
        devtools::document()
        system(paste(
          "R CMD INSTALL --preclean --no-multiarch --with-keep.source",
          getwd()
        ))
      }
    }
  }

  return(sh5HT)

}
