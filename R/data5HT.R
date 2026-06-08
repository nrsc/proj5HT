#' Assemble the dataframes for use
#'
#' @returns dfs
#' @export data5HT
#'
#' @examples
data5HT = function(){

  dfs = list()


  ## Data input from sheets to select for inclusion
  sMdta = readODS::read_ods("den/serotonin/Data/selected_merged_metadata.ods", sheet = "Macaque")
  sMdta = sMdta[!is.na(sMdta$cell_name), ]
  sHdta = readODS::read_ods("den/serotonin/Data/selected_merged_metadata.ods", sheet = "Human")
  sHdta = sHdta[!is.na(sHdta$cell_name), ]
  mMD = bind_rows(sMdta,sHdta)

  #pkl_data =

  dfs$ui = list(
    macaque_ui_data = sMdta,
    human_ui_data = sHdta,
    all_ui_data = mMD,
    ui_cells_MD = headBCI::sheets$MD[headBCI::sheets$MD$cell_name %in% mMD$cell_name,]
  )

  # Compiled dataframes from analysis
  dfs$compiled = list(
    compiled_spikeTTL = read.csv("data-raw/compile_srtPuff.csv"),
  )


  #cell_list_L5_corr = read.csv("data-raw/cell_list_L5_corr.csv")### Don't know where this came from
  #cell_list_L5_corr$patched_cell_container %in% headBCI::sheets$MD$patched_cell_container

  ## Map out the annotation data for cristinas
  dfs$map = list(
    L23Cristina = read.csv("data-raw/HCN23_dataset.csv"),
    smryL23HCN = projHCT::sheets$map$mapCellsNHP[projHCT::sheets$map$mapCellsNHP$cell_name %in% L23Cristina$cell_name,],
    annoL23HCN= projHCT::sheets$map$NHP_anno[projHCT::sheets$map$NHP_anno$cell_name_label %in% L23Cristina$cell_name,],
  )

  #


  return(dfs)

}
