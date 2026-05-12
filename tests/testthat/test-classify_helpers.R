test_that("smooth_vec rollmean preserves length and smooths a step", {
  y <- c(rep(0, 10), rep(100, 10))
  ys <- smooth_vec(y, method = "rollmean", k = 3)
  expect_equal(length(ys), length(y))
  # middle of the smoothed result should be intermediate
  expect_true(any(ys > 0 & ys < 100, na.rm = TRUE))
})

test_that("has_persistent_excursion detects sustained crossings", {
  t <- seq(0, 10, by = 0.1)
  y <- rep(100, length(t))
  y[t >= 2 & t <= 7] <- 130    # 5 s excursion above 115

  expect_true(has_persistent_excursion(t, y, c(0, 10), "exc",
                                       thr = 15, min_persist_s = 1))
  expect_false(has_persistent_excursion(t, y, c(0, 10), "inh",
                                        thr = 15, min_persist_s = 1))
  # not long enough at very high min_persist
  expect_false(has_persistent_excursion(t, y, c(0, 10), "exc",
                                        thr = 15, min_persist_s = 10))
})

test_that("get_extreme_in_window finds min and max", {
  t <- seq(0, 10, by = 0.1)
  y <- 100 + 20 * sin(t)
  hi <- get_extreme_in_window(t, y, c(0, 10), "max")
  lo <- get_extreme_in_window(t, y, c(0, 10), "min")
  expect_gt(hi$value, 110)
  expect_lt(lo$value, 90)
})

test_that("classify_magnitude bins as expected", {
  expect_equal(classify_magnitude(NA_real_),  NA_character_)
  expect_equal(classify_magnitude(5),         "none")
  expect_equal(classify_magnitude(15),        "mild")
  expect_equal(classify_magnitude(30),        "moderate")
  expect_equal(classify_magnitude(80),        "strong")
})

test_that("assert_has_cols errors on missing columns", {
  expect_error(assert_has_cols(data.frame(a = 1), c("a", "b")), "b")
  expect_true(assert_has_cols(data.frame(a = 1, b = 2), c("a", "b")))
})
