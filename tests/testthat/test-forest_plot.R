# tests/testthat/test-forest_plot.R
context("forest_plot")

test_that("forest_plot returns a ggplot object for model summary data", {
  # Create a small dummy summary data frame (as if returned by build_model)
  res_df <- data.frame(
    table_id = c("chr12_pos11161", "chr12_pos12806"),
    outcome  = "phenotype",
    model    = "logistic",
    beta     = c(1.5, 0.0),
    se       = c(0.7, 0.5),
    p        = c(0.05, 1.00),
    OR       = exp(c(1.5, 0.0)),
    OR_lo    = exp(c(1.5 - 1.96*0.7, 0.0 - 1.96*0.5)),
    OR_hi    = exp(c(1.5 + 1.96*0.7, 0.0 + 1.96*0.5))
  )
  # Generate forest plot (use table_id as label, facet by outcome)
  gp <- forest_plot(res_df, label_col = "table_id", facet_col = "outcome", show_errors = TRUE)
  
  # The output should be a ggplot object
  expect_s3_class(gp, "ggplot")
  # The plot should have a vertical reference line and points (check at least one layer for points)
  expect_true(length(gp$layers) >= 1)
})
