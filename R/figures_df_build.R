#' Build the family of dataframes used to render proj5HT figures
#'
#' Loads the compiled puff dataset, joins in cluster metadata, derives a
#' time-binned `x_bin` axis and a unified `assigned_depth`, then carves the
#' master dataframe into the experiment-specific subsets that the figure
#' scripts consume. Each subset is a filtered view of the master `df`; columns
#' are identical across subsets so the same plotting helpers can be reused.
#'
#' @return A named list of data frames. Always contains:
#'   `df`, `df0`, `df_inh`, `dfH`, `df46`, `df_5ct`, `df_way`, `df_ket`,
#'   `df_kcon`, `df_kwo`, `df_carb`, `df_control`, `df_fetal`.
#' @export
figure_df_build <- function() {

  ## ------------------------------------------------------------------ helpers
  # Within each cell, ensure every 2 s bin between the cell's first and last
  # bin is represented. Missing bins correspond to periods of no spiking and
  # are filled with `percent_change = 0` so they participate in the average
  # rather than being silently dropped.
  fill_missing_bins <- function(df, bin_step = 2) {
    df %>%
      dplyr::group_by(cell_name, assigned_subclass, Species) %>%
      tidyr::complete(
        x_bin = seq(min(x_bin), max(x_bin), by = bin_step),
        fill = list(percent_change = 0)
      ) %>%
      dplyr::ungroup()
  }

  ## --------------------------------------------------------------------- data
  MD      <- projHCT::sheets$MD
  comp5HT <- readRDS("data-raw/compiled5HT.rds")

  df <- data.frame(comp5HT$srtPuff)

  ## ----------------------------------------------------------- manual excludes
  bad_cells <- c("Q21.26.014.1A.21.02", "Q21.26.010.11.13.03")
  df <- dplyr::filter(df, !cell_name %in% bad_cells)

  ## ------------------------------------------------------------- normalisation
  df$assigned_subclass <- gsub("/", "", df$assigned_subclass)
  df$Date              <- as.Date(df$Date)

  ## -------------------------------------------------------- derived variables
  # unified depth (LIMs preferred, otherwise from-pia fluorescence proxy)
  df$assigned_depth <- ifelse(is.na(df$LIMs_depth),
                              df$DFP * 0.77,
                              df$LIMs_depth)

  # superficial / deep split at 500 µm from pia
  df$depth_bin <- ifelse(is.na(df$assigned_depth), NA_character_,
                  ifelse(df$assigned_depth < 500, "superficial", "deep"))
  df$depth_bin <- factor(df$depth_bin, levels = c("superficial", "deep"))

  # bring in cluster_Corr from the master metadata sheet
  df <- df %>%
    dplyr::left_join(
      MD %>%
        dplyr::filter(cell_name %in% unique(df$cell_name)) %>%
        dplyr::distinct(cell_name, cluster_Corr),
      by = "cell_name"
    )

  # 2 s time bins for averaging / plotting
  df <- df %>% dplyr::mutate(x_bin = floor(time / 2) * 2)

  df <- fill_missing_bins(df)

  ## -------------------------------------------------------- subclass groupings
  pyr_subclasses  <- c("L23_IT", "L5_ET", "L5_IT",
                       "L6_IT",  "L56_NP", "L6_IT_Car3", "L6_CT", "L6b")
  l46_subclasses  <- c("L6_IT", "L56_NP", "L6_IT_Car3", "L6_CT", "L6b")
  l235_subclasses <- c("L23_IT", "L5_ET", "L5_IT")

  df_inh <- df %>% dplyr::filter(!assigned_subclass %in% pyr_subclasses)
  df     <- df %>% dplyr::filter( assigned_subclass %in% pyr_subclasses)
  df46   <- df %>% dplyr::filter( assigned_subclass %in% l46_subclasses)
  dfH    <- df %>% dplyr::filter(Species == "Human")

  ## ----------------------------------------------------- experiment-level cuts
  ket_contam_window <- c(as.Date("2025-09-17"), as.Date("2026-01-14"))

  df_way     <- df %>% dplyr::filter(bath == "WAY635")
  df_ket     <- df %>% dplyr::filter(bath == "Ket")
  df_5ct     <- df %>% dplyr::filter(puff %in% c("5CT[200nM]", "5CT[50nM]"))
  df_carb    <- df %>% dplyr::filter(puff == "Carb[50uM]")
  df_control <- df %>% dplyr::filter(puff %in% c("Control", "acsf"))
  df_fetal   <- df %>% dplyr::filter(expCon == "Fetal")

  # cells run during the ketanserin contamination window
  df_kcon <- df %>%
    dplyr::filter(Date >= ket_contam_window[1],
                  Date <= ket_contam_window[2])

  # the ketanserin washout cohort
  kwo_cells <- c(
    "QN26.26.004.20.03.01",
    "QN26.26.004.20.03.02",
    "QN26.26.004.20.03.03",
    "QN26.26.004.20.03.04",
    "QN26.26.004.20.05.01",
    "QN26.26.004.20.08.01",
    "QN26.26.004.20.08.02"
  )
  df_kwo <- df %>% dplyr::filter(cell_name %in% kwo_cells)

  # standard puff experiments (no bath drug, outside Ket contam window)
  df0 <- df %>%
    dplyr::filter(
      assigned_subclass %in% l235_subclasses,
      bath   == "none",
      expCon %in% c("Standard_Puff", "lowEGTA_IC"),
      !dplyr::between(Date, ket_contam_window[1], ket_contam_window[2]),
      grepl("^5HT", puff)
    )

  list(
    df         = df,
    df0        = df0,
    df_inh     = df_inh,
    dfH        = dfH,
    df46       = df46,
    df_5ct     = df_5ct,
    df_way     = df_way,
    df_ket     = df_ket,
    df_kcon    = df_kcon,
    df_kwo     = df_kwo,
    df_carb    = df_carb,
    df_control = df_control,
    df_fetal   = df_fetal
  )
}
