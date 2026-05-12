# These tests exercise the full pipeline if the on-disk data fixtures are
# available. They are skipped otherwise so the suite remains portable.

skip_unless_data <- function() {
  rds <- file.path("..", "..", "data-raw", "compiled5HT.rds")
  if (!file.exists(rds)) {
    skip("compiled5HT.rds not available")
  }
  if (!requireNamespace("projHCT", quietly = TRUE)) {
    skip("projHCT not installed")
  }
}

test_that("figure_df_build returns the expected named dataframes", {
  skip_unless_data()

  owd <- setwd(file.path("..", ".."))
  on.exit(setwd(owd), add = TRUE)

  source("R/figures_df_build.R", local = TRUE)
  dfs <- figure_df_build()

  required <- c("df", "df0", "df_inh", "dfH", "df46", "df_5ct",
                "df_way", "df_ket", "df_kcon", "df_kwo", "df_carb",
                "df_control", "df_fetal")
  expect_true(all(required %in% names(dfs)))

  # master df should expose the bin axis and the derived depth columns
  expect_true(all(c("x_bin", "assigned_depth", "depth_bin",
                    "percent_change", "cell_name",
                    "assigned_subclass") %in% names(dfs$df)))
  expect_s3_class(dfs$df$depth_bin, "factor")
  expect_setequal(levels(dfs$df$depth_bin), c("superficial", "deep"))
})

test_that("classify_response5HT runs on the real df0 cohort", {
  skip_unless_data()

  owd <- setwd(file.path("..", ".."))
  on.exit(setwd(owd), add = TRUE)

  source("R/figures_df_build.R", local = TRUE)
  dfs <- figure_df_build()

  per_cell <- classify_response5HT(
    dplyr::filter(dfs$df0,
                  assigned_subclass %in% c("L23_IT", "L5_IT", "L5_ET")),
    rapid_window    = c(0, 10),
    extended_window = c(10, 50),
    thr_pct         = 15,
    min_persist_s   = 1
  )

  expect_gt(nrow(per_cell), 0)
  expect_true("response_type" %in% names(per_cell))
  expect_true(all(as.character(per_cell$response_type) %in%
                    c("excitation", "inhibition",
                      "biphasic_exc_inh", "biphasic_inh_exc",
                      "no_change")))
})
