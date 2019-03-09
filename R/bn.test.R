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
                    alpha = 0.10,
                    longlat = FALSE, noc = TRUE,
                    modified = FALSE) {
  # sanity checking
  arg_check_bn_test(coords, cases, pop, cstar, longlat,
                    alpha, noc)

  coords = as.matrix(coords)
  # estimate of constant risk
  r = sum(cases) / sum(pop)
  # intercentroid distances
  d = sp::spDists(coords, longlat = longlat)

  # find smallest windows with at least c* cases
  cwins = casewin(d, cases, cstar)
  # determine size of each window
  l = sapply(cwins, length)
  # determine population in each window
  case_cwins = sapply(cwins, function(x) sum(cases[x]))
  ncwins = sapply(cwins, function(x) sum(pop[x]))
  ex = r * ncwins

  if (!modified) {
    pvalue = stats::ppois(cstar - 1, lambda = ex,
                          lower.tail = FALSE)
    } else {
    pvalue = stats::ppois(case_cwins - 1, lambda = ex,
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
            pop = pop, ex = r * pop, longlat = longlat,
            w = NULL, d = d)
}

#' Argument checking for bn.test
#'
#' Check the arguments of the bn.test function
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
