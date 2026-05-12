test_that("classify_response5HT calls each canonical shape correctly", {
  df <- make_synth_df()

  per_cell <- classify_response5HT(
    df,
    rapid_window    = c(0, 10),
    extended_window = c(10, 30),
    thr_pct         = 15,
    min_persist_s   = 1,
    smooth_k        = 3
  )

  expect_s3_class(per_cell, "tbl_df")
  expect_setequal(per_cell$cell_name,
                  c("c_flat", "c_exc", "c_inh", "c_ei", "c_ie"))

  call_for <- function(cn) {
    as.character(per_cell$response_type[per_cell$cell_name == cn])
  }

  expect_equal(call_for("c_flat"), "no_change")
  expect_equal(call_for("c_exc"),  "excitation")
  expect_equal(call_for("c_inh"),  "inhibition")
  expect_equal(call_for("c_ei"),   "biphasic_exc_inh")
  expect_equal(call_for("c_ie"),   "biphasic_inh_exc")
})

test_that("classify_response5HT carries metadata and exposes window flags", {
  df <- make_synth_df()
  per_cell <- classify_response5HT(df, smooth_k = 3, min_persist_s = 1)

  expect_true(all(c("rapid_exc", "rapid_inh", "ext_exc", "ext_inh",
                    "exc_peak_time", "inh_peak_time") %in% names(per_cell)))
  expect_true("assigned_subclass" %in% names(per_cell))
  expect_true("Species" %in% names(per_cell))

  exc_row <- per_cell[per_cell$cell_name == "c_exc", ]
  expect_true(exc_row$rapid_exc)
  expect_false(exc_row$rapid_inh)

  ei_row <- per_cell[per_cell$cell_name == "c_ei", ]
  expect_true(ei_row$rapid_exc)
  expect_true(ei_row$ext_inh)
})

test_that("classify_response5HT errors on missing required columns", {
  bad <- data.frame(cell_name = "x", time = 1:5)
  expect_error(classify_response5HT(bad), "percent_change")
})
