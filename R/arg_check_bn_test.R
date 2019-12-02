#' Argument checking for bn.test
#'
#' Check the arguments of the bn.test function
#' @param coords An \eqn{n \times 2} matrix of centroid
#'   coordinates for the regions.
#' @param cases The number of cases observed in each region.
#' @param pop The population size associated with each
#'   region.
#' @param nsim The number of simulations from which to
#'   compute the p-value.
#' @param alpha The significance level to determine whether
#'   a cluster is signficant.  Default is 0.10.
#' @param longlat The default is \code{FALSE}, which
#'   specifies that Euclidean distance should be used.If
#'   \code{longlat} is \code{TRUE}, then the great circle
#'   distance is used to calculate the intercentroid
#'   distance.
#' @param noc A logical value indicating whether all
#'   significant clusters should be returned (\code{FALSE})
#'   or only the non-overlapping clusters (\code{TRUE})
#'   arranged in order of significance.  The default is
#'   \code{TRUE}.
#' @param modified A logical value indicating whether a
#'   modified version of the test should be performed.  The
#'   original paper recommends computing the p-value for
#'   each cluster as \code{1 - ppois(cstar - 1, lambda =
#'   expected)}. The modified version replaces \code{cstar}
#'   with \code{cases}, the observed number of cases in the
#'   region, and computes the p-value for the cluster as
#'   \code{1 - ppois(cases - 1, lambda = ex)}. The default
#'   is \code{modified = FALSE}.
#' @param cstar A non-negative integer indicating the
#'   minimum number of cases to include in each window.
#' @return NULL
#' @export
#' @keywords internal
arg_check_bn_test = function(coords, cases, pop, cstar,
                             longlat, alpha, noc) {
  if (!(is.matrix(coords) | is.data.frame(coords))) {
    stop("coords should be a matrix or a data frame")
  }
  if (ncol(coords) != 2) {
    stop("coords must have two columns")
  }
  N = nrow(coords)
  if (length(cases) != N) {
    stop("length(cases) != nrow(coords)")
  }
  if (!is.numeric(cases)) {
    stop("cases should be a numeric vector")
  }
  if (length(pop) != N) {
    stop("length(pop) != nrow(coords)")
  }
  if (!is.numeric(pop)) {
    stop("pop should be a numeric vector")
  }
  if (length(cstar) != 1 || !is.numeric(cstar)) {
    stop("cstar should be a numeric vector of length 1")
  }
  if (cstar < 1 || cstar > sum(cases)) {
    stop("cstar should be at least 1 and less than or equal to the sum(cases)")
  }
  if (length(longlat) != 1) {
    stop("length(longlat) != 1")
  }
  if (!is.logical(longlat)) {
    stop("longlat should be a logical value")
  }
  if (length(alpha) != 1 || !is.numeric(alpha)) {
    stop("alpha should be a numeric vector of length 1")
  }
  if (alpha < 0 || alpha > 1) {
    stop("alpha should be a value between 0 and 1")
  }
  if (length(noc) != 1) {
    stop("length(noc) != 1")
  }
  if (!is.logical(noc)) {
    stop("noc should be a logical value")
  }
}
