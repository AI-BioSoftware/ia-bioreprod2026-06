test_that("prepare_figure2a_data returns expected structure", {
  prepared <- prepare_figure2a_data(default_excel_path())
  expect_type(prepared, "list")
  expect_named(prepared, c("data", "expr_z", "time_vals", "gene_ids"))
  expect_equal(nrow(prepared$expr_z), 1246)
  expect_equal(length(prepared$time_vals), ncol(prepared$expr_z))
  expect_true(all(diff(prepared$time_vals) >= 0))
})

test_that("selected genes all satisfy LS cutoff", {
  prepared <- prepare_figure2a_data(default_excel_path())
  expect_true(all(prepared$data$LS_cutoff == "Yes"))
})

test_that("z-score matrix is finite", {
  prepared <- prepare_figure2a_data(default_excel_path())
  expect_true(all(is.finite(prepared$expr_z)))
})

test_that("missing file triggers a clear error", {
  expect_error(
    prepare_figure2a_data("does-not-exist.xlsx"),
    "Input file does not exist"
  )
})
