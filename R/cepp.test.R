#' Cluster Evalation Permutation Procedure Test
#'
#' \code{cepp.test} implements the Cluster Evaluation
#' Permutation Procedure test of Turnbull et al. (1990)
#' for finding disease clusters.
#'
#' @inheritParams scan.test
#' @inheritParams casewin
#' @inheritParams bn.test
#' @param nstar The size of the at-risk population
#' in each window.
#' @param simdist A character string indicating whether the
#' simulated data should come from a \code{"multinomial"}
#' or \code{"poisson"} distribution.  The default is
#' \code{"multinomial"}, which fixes the total number of
#' cases observed in each simulated data set.
#'
#' @return Returns a \code{scan} object.
#' @author Joshua French
#' @seealso \code{\link{scan.stat}},
#'   \code{\link{plot.scan}}, \code{\link{scan.test}}
#' @export
#' @references Bruce W. Turnbull, Eric J. Iwano, William S. Burnett,
#'   Holly L. Howe, Larry C. Clark (1990).  Monitoring for Clusters of Disease:
#'   Application to Leukemia Incidence in Upstate New York,
#'   American Journal of Epidemiology, 132(supp1):136-143.
#'   <doi:10.1093/oxfordjournals.aje.a115775>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(x, y))
#' cases = nydf$cases
#' pop = nydf$pop
#' out = cepp.test(coords = coords, cases = cases, pop = pop,
#'                 nstar = 1000, alpha = 0.1)
#' plot(out)
#'
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
cepp.test = function(coords, cases, pop, nstar,
                     ex = sum(cases) / sum(pop) * pop,
                     nsim = 499, alpha = 0.10,
                     longlat = FALSE, noc = TRUE,
                     simdist = "multinomial") {
  # sanity checking
  arg_check_cepp_test(coords = coords, cases = cases,
                      pop = pop, nstar = nstar, ex = ex,
                      nsim = nsim, alpha = alpha,
                      longlat = longlat, noc = noc,
                      simdist = simdist)

  coords = as.matrix(coords)

  # intercentroid distances
  d = sp::spDists(coords, longlat = longlat)

  # find smallest windows with at least n* pop
  nn = casewin(d, pop, nstar)

  # determine nn weights
  wts = cepp.weights(nn, pop, nstar)

  # observed number of cases in each window
  cstar = sapply(seq_along(nn), function(i) {
    sum(cases[nn[[i]]] * wts[[i]])
  }, USE.NAMES = FALSE)

  csim = cepp.sim(nsim = nsim, nn = nn, ty = sum(cases),
                  ex = ex, wts = wts, simdist = simdist)

  pvalue = mc.pvalue(cstar, csim)

  op = order(cstar, decreasing = TRUE)

  if (noc) {
    # determine idx of unique non-overlapping clusters in
    # order of significance
    u = smacpod::noc(nn[op])
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
  sig_regions = nn[op]
  sig_tstat = cstar[op]
  sig_p = pvalue[op]
  prep.scan2(tobs = sig_tstat, zones = sig_regions,
            pvalue = sig_p, coords = coords, cases = cases,
            pop = pop, ex = ex, longlat = longlat,
            w = NULL, d = d)
}
