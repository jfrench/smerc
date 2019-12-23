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
#' @param ex A vector of expected counts
#' @param modified A logical value for whether the test should be "modified"
#' @noRd
arg_check_bn_test = function(coords, cases, pop, cstar,
                             longlat, alpha, ex,
                             modified) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_cstar(cstar, cases)
  arg_check_longlat(longlat)
  arg_check_alpha(alpha)
  arg_check_ex(ex, N)
  arg_check_modified(modified)
}
