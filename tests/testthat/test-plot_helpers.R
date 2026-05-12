test_that("plot_puff_traces returns a ggplot with grey + black layers", {
  source(file.path("..", "..", "figs", "scripts", "plot_helpers.R"),
         local = TRUE)
  df <- make_synth_df()

  p <- plot_puff_traces(df,
                        xlim = c(-10, 30),
                        ylim = c(0, 250),
                        bin_step = 2,
                        min_cells_per_bin = 1)

  expect_s3_class(p, "ggplot")

  # two geom_line layers: grey per-cell traces + black mean
  line_layers <- vapply(p$layers,
                        function(l) inherits(l$geom, "GeomLine"),
                        logical(1))
  expect_gte(sum(line_layers), 2)

  # rendering should not error
  expect_silent(ggplot2::ggplot_build(p))
})

test_that("plot_puff_traces_by_response builds when classifier is run", {
  source(file.path("..", "..", "figs", "scripts", "plot_helpers.R"),
         local = TRUE)
  df <- make_synth_df()
  per_cell <- classify_response5HT(df, smooth_k = 3, min_persist_s = 1)

  p <- plot_puff_traces_by_response(df, per_cell,
                                    xlim = c(-10, 30),
                                    ylim = c(0, 250),
                                    bin_step = 2,
                                    min_cells_per_bin = 1)
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(p))
})
