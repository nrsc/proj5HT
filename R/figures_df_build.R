figure_df_build = function(){


  fill_missing_bins <- function(df,
                                bin_col = "x_bin",
                                value_col = "percent_change",
                                bin_seq = seq(-10, 50, by = 2)) {

    df %>%
      dplyr::group_by(cell_name, assigned_subclass, Species) %>%
      tidyr::complete(
        !!rlang::sym(bin_col) := bin_seq,
        fill = setNames(list(0), value_col)
      ) %>%
      dplyr::ungroup()
  }

  #comp5HT = compile5HT()
  MD = projHCT::sheets$MD
  comp5HT = readRDS("data-raw/compiled5HT.rds")

  #### Master Dataframe prep ####
  df = data.frame(comp5HT$srtPuff)

  #### Manual filter of experiments ####
  df = df %>%  dplyr::filter(
      !cell_name %in% c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03")
    )


  #### Manipulate already present variables if necessary ####
  df$expCon = gsub("standard_puff", "Standard_Puff", df$expCon)
  df$assigned_subclass = gsub("/", "", df$assigned_subclass)
  df$Date = as.Date(df$Date)

  #### Add variables ####
  ### Depth ##
  df$assigned_depth = ifelse(is.na(df$LIMs_depth), df$DFP*.77, df$LIMs_depth)
  ### Cluster_corr binding ###
  df <- df %>%
    dplyr::left_join(MD[MD$cell_name %in% unique(df$cell_name), c("cluster_Corr", "cell_name")] %>%
                       dplyr::distinct(cell_name, cluster_Corr), by = "cell_name")
  #### Create 2 second time bins ####
  df = df %>% mutate(x_bin = floor(time / 2) * 2)

  # #### Fill in empty values ####
  # df <- df %>%
  #   dplyr::group_by(cell_name,
  #                   assigned_subclass,
  #                   Species) %>%
  #   tidyr::complete(
  #     x_bin = seq(min(x_bin), max(x_bin), by = 2),
  #     fill = list(percent_change = 0)
  #   ) %>%
  #   dplyr::ungroup()
  #
  # df %>%
  #   dplyr::group_by(x_bin, cell_name, assigned_subclass, Species) %>%
  #   dplyr::summarise(
  #     y = mean(percent_change),
  #     .groups = "drop"
  #   )

  df <- fill_missing_bins(df)

  #### Subclass specific dataframes ####
  df = df %>% dplyr::filter(
    assigned_subclass %in% c("L6_IT", "L56_NP", "L6_IT_Car3","L6_CT", "L6b","L23_IT", "L5_ET", "L5_IT")
  )

  dfH <- df[df$Species == "Human", ]

  #### Overview of conditions ####
  # table(df$bath)
  # table(df$expCon)
  # table(df$puff)

  #### Bath drug data from masterDF ####
  df_way = df %>%
    dplyr::filter(
      bath == "WAY635")

  df_ket = df %>%
    dplyr::filter(
      bath == "Ket")

  #### Isolate sparse datasets
  #Ket conamination period
  df_kcon = df %>%
    dplyr::filter(
      Date >= as.Date("2025-09-17") &
        Date <= as.Date("2026-01-14")
    )

  #Ket washout experiment
  df_kwo = df %>%
    dplyr::filter(
      cell_name %in% c(
        "QN26.26.004.20.03.01",
        "QN26.26.004.20.03.02",
        "QN26.26.004.20.03.03",
        "QN26.26.004.20.03.04",
        "QN26.26.004.20.05.01",
        "QN26.26.004.20.08.01",
        "QN26.26.004.20.08.02"
      )
    )


  #### L46 Subclass specific dataframes ####
  df46 = df %>% dplyr::filter(
    assigned_subclass %in% c("L6_IT", "L56_NP", "L6_IT_Car3","L6_CT", "L6b")
  )


  #### Get standard puff experiments ####
  df0 = df %>%
    dplyr::filter(
      assigned_subclass %in% c("L23_IT", "L5_ET", "L5_IT"),
      bath == "none",
      expCon %in% c("Standard_Puff", "lowEGTA_IC"),
      !dplyr::between(
        Date,
        as.Date("2025-09-17"),
        as.Date("2026-01-14")
      ),
      grepl("^5HT", puff))

  dfs = list(
    df = df,
    df0 = df0,
    dfH = dfH,
    df46 = df46,
    df_way = df_way,
    df_ket = df_ket,
    df_kcon = df_kcon
  )

  return(dfs)

}
