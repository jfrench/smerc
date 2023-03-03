if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(smerc)
  test_check("smerc")
} else {
  message("tests could not be run because testthat is not installed")
}
