if (requireNamespace("lintr", quietly = TRUE)) {
  context("lints")
  mylints = lintr::with_defaults(assignment_linter = NULL,
                                 closed_curly_linter = NULL,
                                 commas_linter = NULL)
  test_that("Package Style", {
    lintr::expect_lint_free(linters = mylints)
  })
}
