#' Compute middle p-value
#'
#' Computes P(Y > cases) + P(Y = cases)/2 when Y ~
#' Poisson(ex).  This is middle p-value computed by Tango
#' and Takahashi (2012)
#'
#' @inheritParams scan.test
#' @return A vector of middle p-values
#' @export
#' @author Joshua French
#' @references Tango, T. and Takahashi, K. (2012), A
#'   flexible spatial scan statistic with a restricted
#'   likelihood ratio for detecting disease clusters.
#'   Statist. Med., 31: 4207-4218. <doi:10.1002/sim.5478>
#'
#' @examples
#' rflex.midp(5, 3)
#' 1 - ppois(5, 3) + dpois(5, 3)/2
rflex.midp = function(cases, ex) {
  if (length(cases) != length(ex)) {
    stop("length(cases) != length(ex)")
  }
  if (!is.numeric(cases) | !is.numeric(ex)) {
    stop("cases and ex must be numeric vectors")
  }
  stats::ppois(cases, ex, lower.tail = FALSE) + 
    stats::dpois(cases, ex)/2
}