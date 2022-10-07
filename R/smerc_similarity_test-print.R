#' Print object of class \code{smerc_similarity_test}.
#'
#' Print a \code{smerc_similarity_test} object. If the \code{crayon} package
#' is installed, then the results are printed in color.
#'
#' @param x An object of class \code{smerc_similarity+test}.
#' @param digits Number of significant digits to print.
#' @param \dots Not currently implemented.
#' @export
print.smerc_similarity_test <- function(x, ..., digits = 2) {
  cat("\n")
  cat("Walter (1992) Constant-risk Moran's I\n")
  cat("\n")
  cat("null hypothesis: rho = 0", "\n")
  if (x$alternative == "two.sided") {
    cat("alternative hypothesis: rho != 0", "\n")
  } else if (x$alternative == "greater") {
    cat("alternative hypothesis: rho > 0", "\n")
  } else if (x$alternative == "less") {
    cat("alternative hypothesis: rho < 0", "\n")
  }
  cat("test statistic:", round(x$test_statistic, digits = digits), "\n")
  cat("p-value:", round(x$pvalue, digits = digits), "\n")
  cat("nsim:", x$nsim, "\n")
  cat("simulation procedure: ", x$simdist, "\n")
}
