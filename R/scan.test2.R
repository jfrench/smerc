#' Spatial Scan Test
#'
#' \code{scan.test} performs the original spatial scan test
#' of Kulldorf (1997) based on a fixed number of cases.
#' Candidate zones are circular and extend from the observed
#' region centroids.  The clusters returned are
#' non-overlapping, ordered from most significant to least
#' significant.  The first cluster is the most  likely to be
#' a cluster.  If no significant clusters are found, then
#' the most likely cluster is returned (along with a
#' warning).
#'
#' @param coords An \eqn{n \times 2} matrix of centroid
#'   coordinates for the regions in the form (x, y) or
#'   (longitude, latitude) is using great circle distance.
#' @param cases The number of cases observed in each region.
#' @param pop The population size associated with each
#'   region.
#' @param ex The expected number of cases for each region.
#'   The default is calculated under the constant risk
#'   hypothesis.
#' @param nsim The number of simulations from which to
#'   compute the p-value.
#' @param ubpop The upperbound of the proportion of the
#'   total population to consider for a cluster.
#' @param alpha The significance level to determine whether
#'   a cluster is signficant.  Default is 0.10.
#' @param longlat The default is \code{FALSE}, which
#'   specifies that Euclidean distance should be used. If
#'   \code{longlat} is \code{TRUE}, then the great circle
#'   distance is used to calculate the intercentroid
#'   distance.
#' @param type The type of scan statistic to compute. The
#'   default is \code{"poisson"}. The other choice
#'   is \code{"binomial"}.
#' @param min.cases The minimum number of cases required for
#'   a cluster.  The default is 2.
#' @param simdist Character string indicating the simulation
#' distribution. The default is \code{"multinomial"}, which
#' conditions on the total number of cases observed. The
#' other options are \code{"poisson"} and \code{"binomial"}
#' @inheritParams pbapply::pblapply
#'
#' @return Returns a \code{smerc_cluster} object.
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial scan
#'   statistic. Communications in Statistics - Theory and
#'   Methods, 26(6): 1481-1496,
#'   <doi:10.1080/03610929708831995>
#'
#' Waller, L.A. and Gotway, C.A. (2005). Applied Spatial
#' Statistics for Public Health Data. Hoboken, NJ: Wiley.
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 499,
#'                 alpha = 0.10, longlat = TRUE)
#' out2 = scan.test2(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 499,
#'                 alpha = 0.10, longlat = TRUE)
scan.test2 = function(coords, cases, pop,
                     ex = sum(cases) / sum(pop) * pop,
                     nsim = 499, alpha = 0.1,
                     ubpop = 0.5, longlat = FALSE, cl = NULL,
                     type = "poisson",
                     min.cases = 2,
                     simdist = "multinomial") {
  # argument checking
  type = match.arg(type, c("poisson", "binomial"))
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_scan_test(coords = coords, cases = cases,
                      pop = pop, ex = ex, nsim = nsim,
                      alpha = alpha, ubpop = ubpop,
                      longlat = longlat,
                      k = 1, w = diag(nrow(coords)),
                      type = type, simdist = simdist,
                      min.cases = min.cases)

  # convert to proper format
  coords = as.matrix(coords)
  N = nrow(coords)
  # compute inter-centroid distances
  d = sp::spDists(coords, longlat = longlat)

  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = scan.nn(d, pop, ubpop)
  nnn = unlist(lapply(nn, length), use.names = FALSE)

  # determine total number of cases in each successive
  # window, total number of cases
  yin = nn.cumsum(nn, cases)
  ty = sum(cases) # sum of all cases

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = nn.cumsum(nn, ex)
    eout = sum(ex) - ein
    # correct for the situation when the expected number of cases
    # is not the same as the observed number of cases
    mult = ty / sum(ex)
    ein = ein * mult
    eout = eout * mult
    popin = NULL
    popout = NULL
    tpop = NULL
    tobs = stat_poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    tpop = sum(pop)
    popin = nn.cumsum(nn, pop)
    popout = tpop - popin
    tobs = stat_binom(yin, ty - yin, ty, popin, popout, tpop)
  }
  # tobs in nn format
  tobs_nn = split(tobs, f = rep(seq_along(nn), times = nnn))

  noc_info = noc_nn(nn, tobs_nn)
  tobs = noc_info$tobs

  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = scan.sim(nsim = nsim, nn = nn, ty = ty,
                    ex = ex, type = type, ein = ein,
                    eout = eout, popin = popin,
                    popout = popout, tpop = tpop, cl = cl,
                    simdist = simdist, pop = pop, min.cases = min.cases)
    pvalue = mc.pvalue(tobs, tsim)
  } else {
    pvalue = rep(1, length(tobs))
  }

  # significant, ordered, non-overlapping clusters and
  # information
  pruned = sig_prune(tobs = tobs, zones = noc_info$clusts,
                     pvalue = pvalue, alpha = alpha)

  smerc_cluster(tobs = pruned$tobs, zones = pruned$zones,
                pvalue = pruned$pvalue, coords = coords,
                cases = cases, pop = pop, ex = ex,
                longlat = longlat, method = "circular scan",
                rel_param = list(type = type,
                                 simdist = simdist,
                                 nsim = nsim,
                                 ubpop = ubpop,
                                 min.cases = min.cases),
                alpha = alpha,
                w = NULL, d = d)
}
