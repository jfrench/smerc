#' Besag-Newell Test
#'
#' \code{bn.test} implements the Besag-Newell test of Besag
#' and Newell (1991) for finding disease clusters.
#'
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
#' @inheritParams scan.test
#' @inheritParams casewin
#'
#' @return Returns a \code{scan} object.
#' @author Joshua French
#' @seealso \code{\link{scan.stat}},
#'   \code{\link{plot.scan}}, \code{\link{scan.test}},
#'   \code{\link{flex.test}}, \code{\link{dmst.test}},
#'   \code{\link{uls.test}}, \code{\link{mlf.test}}
#' @export
#' @references Besag, J. and Newell, J.  (1991). The
#'   detection of clusters in rare diseases, Journal of the
#'   Royal Statistical Society, Series A, 154, 327-333.
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(x, y))
#' out = bn.test(coords = coords, cases = nydf$cases,
#'               pop = nydf$pop, cstar = 6,
#'               alpha = 0.1)
#' plot(out)
#'
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
bn.test = function(coords, cases, pop, cstar,
                   ex = sum(cases) / sum(pop) * pop,
                   alpha = 0.10,
                   longlat = FALSE, noc = TRUE,
                   modified = FALSE) {
  # sanity checking
  arg_check_bn_test(coords, cases, pop, cstar, longlat,
                    alpha, noc, ex)

  coords = as.matrix(coords)
  # estimate of constant risk
  # r = sum(cases) / sum(pop)
  # intercentroid distances
  d = sp::spDists(coords, longlat = longlat)

  # find smallest windows with at least c* cases
  cwins = casewin(d, cases, cstar)
  # determine size of each window
  l = sapply(cwins, length)
  # determine cases and expected in each window
  case_cwins = zones.sum(cwins, cases)
  ex_cwins = zones.sum(cwins, ex)

  if (!modified) {
    pvalue = stats::ppois(cstar - 1, lambda = ex_cwins,
                          lower.tail = FALSE)
    } else {
    pvalue = stats::ppois(case_cwins - 1, lambda = ex_cwins,
                          lower.tail = FALSE)
  }

  op = order(pvalue)

  if (noc) {
    # determine idx of unique non-overlapping clusters in
    # order of significance
    u = smacpod::noc(cwins[op])
    op = op[u]
    # return only significant clusters
    if (pvalue[op][1] > alpha) {
      warning("No significant clusters.  Returning most likely cluster.")
      op = op[1]
    } else {
      op = op[which(pvalue[op] <= alpha)]
    }
  }

  # for the unique, non-overlapping clusters in order of
  # significance, find the associated test statistic,
  # p-value, centroid, window radius, cases in window,
  # expected cases in window, population in window,
  # standarized mortality ratio, relative risk
  sig_regions = cwins[op]
  sig_tstat = l[op]
  sig_p = pvalue[op]
  prep.scan2(tobs = sig_tstat, zones = sig_regions,
            pvalue = sig_p, coords = coords, cases = cases,
            pop = pop, ex = ex, longlat = longlat,
            w = NULL, d = d)
}
