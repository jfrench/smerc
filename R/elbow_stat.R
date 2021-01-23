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
#'   coordinates for the regions.
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
#'                 pop = nydf$pop, nsim = 0,
#'                 alpha = 1, longlat = TRUE)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "county", region = "new york",
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(countyMapEnv, package = "maps")
#' plot(out, usemap = TRUE, mapargs = mapargs)
#'
#' # a second example to match the results of Waller and Gotway (2005)
#' # in chapter 7 of their book (pp. 220-221).
#' # Note that the 'longitude' and 'latitude' used by them has
#' # been switched.  When giving their input to SatScan, the coords
#' # were given in the order 'longitude' and 'latitude'.
#' # However, the SatScan program takes coordinates in the order
#' # 'latitude' and 'longitude', so the results are slightly different
#' # from the example above.
#' coords = with(nydf, cbind(y, x))
#' out2 = scan.test(coords = coords, cases = floor(nydf$cases),
#'                   pop = nydf$pop, nsim = 0,
#'                   alpha = 1, longlat = TRUE)
#' # the cases observed for the clusters in Waller and Gotway: 117, 47, 44
#' # the second set of results match
#' sget(out2$clusters, name = "cases")[1:3]
elbow_stat = function(coords, cases, pop,
                     ex = sum(cases) / sum(pop) * pop,
                     nsim = 499, alpha = 0.05,
                     ubpop_seq = seq(0.01, 0.5, len = 50),
                     longlat = FALSE, cl = NULL,
                     type = "poisson",
                     min.cases = 2,
                     simdist = "multinomial") {
  # argument checking
  type = match.arg(type, c("poisson", "binomial"))
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  # arg_check_scan_test(coords = coords, cases = cases,
  #                     pop = pop, ex = ex, nsim = nsim,
  #                     alpha = alpha, ubpop = ubpop,
  #                     longlat = longlat,
  #                     k = 1, w = diag(nrow(coords)),
  #                     type = type, simdist = simdist,
  #                     min.cases = min.cases)

  # convert to proper format
  coords = as.matrix(coords)
  N = nrow(coords)
  # compute inter-centroid distances
  d = sp::spDists(coords, longlat = longlat)

  # total number of cases
  ty = sum(cases)
  # total population
  tpop = sum(pop)

  # for each region, determine sorted nearest neighbors
  # subject to LARGEST population constraint
  nn = scan.nn(d, pop, max(ubpop_seq))
  # determine duplicate zones
  wdup = nndup(nn, N)

  # compute number of cases in each candidate zone
  # only keep non-duplicated zones
  yin = nn.cumsum(nn, cases)[!wdup]
  # compute number of cases in each candidate zone
  # only keep non-duplicated zones
  popin = nn.cumsum(nn, pop)[!wdup]
  # determine non-duplicated zones
  zones = nn2zones(nn)[!wdup]

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = nn.cumsum(nn, ex)
    eout = sum(ex) - ein
    popout = NULL
    tobs = stat.poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    popout = tpop - popin
    tobs = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
  }

  # determine tobs for each population upper bound
  tobs_seq = tobs_keep_seq(tobs = tobs,
                           ubpop_seq = ubpop_seq,
                           popin = popin,
                           tpop = tpop,
                           yin = yin,
                           min.cases = min.cases)

  # determine zones for each population upper bound
  zones_seq = tobs_keep_seq(tobs = zones,
                            ubpop_seq = ubpop_seq,
                            popin = popin,
                            tpop = tpop,
                            yin = yin,
                            min.cases = min.cases)

  # compute test statistics for simulated data for all candidate zones
  message("computing statistics for simulated data:")
  tall_yin_sim = scan.sim(nsim = nsim, nn = nn, ty = ty,
                  ex = ex, type = type, ein = ein,
                  eout = eout, popin = popin,
                  popout = popout, tpop = tpop, cl = cl,
                  simdist = simdist, pop = pop,
                  min.cases = 0,
                  return_type = "all")

  # for each population upper bound, determine the maximum test
  # statistic for each simulated data set
  message("partitioning statistics by ubpop_seq:")
  tsim_seq = tsim_ubpop_seq(ubpop_seq = ubpop_seq,
                            tall_yin_sim = tall_yin_sim,
                            min.cases = min.cases,
                            popin = popin,
                            tpop = tpop)
  message("partitioning statistics by ubpop_seq:")
  tsim_seq = tsim_ubpop_seq(ubpop_seq = ubpop_seq,
                            tall_yin_sim = tall_yin_sim,
                            min.cases = min.cases,
                            popin = popin,
                            tpop = tpop)

  message("computing p-values by ubpop_seq")
  pvalue_seq = pbapply::pbmapply(mc.pvalue, tobs = tobs_seq, tsim = tsim_seq)

  # determine significant clusters for each population upper bound
  pruned_seq = suppressWarnings(
    mapply(
      FUN = sig_noc,
      tobs = tobs_seq,
      zones = zones_seq,
      pvalue = pvalue_seq,
      MoreArgs = list(alpha = alpha,
                      order_by = "tobs"),
      SIMPLIFY = FALSE
    )
  )
  # extract the significant (or most likely cluster) tobs
  # for each population upper bound
  sig_tobs_seq = lget(pruned_seq, name = "tobs")

  # compute the negative sum of the test statistics (or MLC)
  neg_lrt = -sapply(sig_tobs_seq, sum)

  # for the significant clusters, compute the total cases and total ex
  # get sig pvalues for each population upper bound
  sig_pvalue_seq = lget(pruned_seq, name = "pvalue")
  # get sig zones for each population upper bound
  sig_zones_seq = lget(pruned_seq, name = "zones")

  # get number of cases in each sig zone for each population upper bound
  sig_yin_seq = lapply(sig_zones_seq, zones.sum, y = cases)
  # get expected cases in each sig zone for each population upper bound
  sig_ein_seq = lapply(sig_zones_seq, zones.sum, y = ex)
  # get order
  sig_order_seq = lapply(sig_ein_seq, order, decreasing = FALSE)

  # # # get sum of cases for all sig zones for each population upper bound
  # sig_cumsum_yin_seq = sapply(sig_yin_seq_zones, cumsum)
  # # get expected cases in each sig zone for each population upper bound
  # sig_ein_seq_zones = sapply(sig_zones_seq, zones.sum, y = ex)
  # # get sum of ex for all sig zones for each population upper bound
  # sig_cumsum_ein_seq = sapply(sig_ein_seq_zones, cumsum)

  # assign any non-significant sums (in case only a MLC returned)
  all_sig = (sapply(pvalue_seq, min) < alpha)


  structure(list(neg_lrt = neg_lrt,
                 ubpop_seq = ubpop_seq),
            class = "smerc_elbow_stats")
}

# for each population upper bound, determine the zones with
# population size <= to that upper bound
#' Title
#'
#' @param popin vector of population in each candidate zone
#' @param tpop total population across all regions
#' @param ubpop_seq # sequence of population upper bounds to consider
#'
#' @return logical matrix indicating which candidate zones
#' satisfy the sequence of population constraints
seq_pop_idx = function(popin, tpop, ubpop_seq) {
  sapply(ubpop_seq, function(ubpop, prop_popin) {
    prop_popin <= ubpop
  }, prop_popin = popin/tpop)
}

# keep the tobs that satisfy two constraints:
# 1. the population size constrained by ubpop
# 2. the number of cases constrained by min.cases
tobs_keep = function(tobs, ubpop, popin, tpop, yin, min.cases) {
  # determine which candidate zones satisfy population constraint
  keep_pop = (popin/tpop <= ubpop)
  # determine which candidate zones satisfy min.cases constraint
  keep_cases = (yin >= min.cases)
  # return the observed statistics that satisfy both requirements
  tobs[keep_pop & keep_cases]
}

# return a list with the tobs that satisfy two constraints:
# 1. the population size constrained by ubpop
# 2. the number of cases constrained by min.cases
# for each ubpop in ubpop_seq
tobs_keep_seq = function(tobs, ubpop_seq, popin, tpop, yin, min.cases) {
  lapply(seq_along(ubpop_seq), function(i) {
    tobs_keep(tobs = tobs,
              ubpop = ubpop_seq[i],
              popin = popin,
              tpop = tpop,
              yin = yin,
              min.cases = min.cases)
  })
}

# based on a specific ubpop, determine the maximum test statistic
# for each simulated data set while accounting for the min.cases
# requirement
tsim_ubpop = function(ubpop, tall_yin_sim, min.cases, popin, tpop) {
  # determine which candidate zones satisfy population constraint
  keep_pop = (popin/tpop <= ubpop)
  # for each simulated data set, determine the maximum test
  # statistic for the zones satisfying the population constraint
  # and the min.cases constraint
  sapply(seq_along(tall_yin_sim), function(i) {
    # determine zones satisfying the min.cases constraint
    keep_cases = (tall_yin_sim[[i]]$yin >= min.cases)
    # max of statistic satisfying both constraints
    max(tall_yin_sim[[i]]$tall[keep_pop & keep_cases])
  })
}

# based on the sequence ubpop_seq, determine the maximum test statistic
# for each simulated data set while accounting for the min.cases
# requirement
tsim_ubpop_seq = function(ubpop_seq, tall_yin_sim, min.cases, popin, tpop) {
  # for each ubpop in ubpop_seq, compute the maximum test statistic
  # for each simulated data set
  pbapply::pbsapply(X = ubpop_seq,
         FUN = tsim_ubpop,
         tall_yin_sim = tall_yin_sim,
         min.cases = min.cases,
         popin = popin,
         tpop = tpop, simplify = FALSE)
}


#' Argument checking for scan tests
#'
#' @param coords A matrix of coordinates
#' @param cases A vector of numeric cases
#' @param pop A vector of population values
#' @param ex A vector of expected counts
#' @param nsim A non-negative integer
#' @param alpha A value greater than 0
#' @param nreport Not used
#' @param ubpop_seq An ascending sequence of values between min(pop)/sum(pop) and 0.5.
#' @param longlat A logical. TRUE is great circle distance.
#' @param parallel Not used.
#' @param k Number of nearest neighbors. Not always needed.
#' @param w A spatial proximity matrix
#' @param type Statistic type
#' @param simdist Distribution of simulation
#' @param min.cases Minimum number of cases. Only for scan.test.
#' @return NULL
#' @noRd
arg_check_elbow_stat =
  function(coords, cases, pop, ex, nsim, alpha,
           nreport = NULL,
           ubpop_seq, longlat, parallel = NULL, k, w, type = NULL,
           simdist = NULL, min.cases = NULL, return_type) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_ex(ex, N)
  arg_check_nsim(nsim)
  arg_check_alpha(alpha)
  # nreport no check, deprecated
  arg_check_ubpop(ubpop)
  arg_check_longlat(longlat)
  # parallel no check, deprecated
  arg_check_k(k, N)
  arg_check_w(w, N)
  if (!is.null(type)) {
    arg_check_type(type)
  }
  if (!is.null(simdist)) {
    arg_check_simdist(simdist)
  }
  arg_check_simdist(simdist)
  if (!is.null(min.cases)) {
    arg_check_min_cases(min.cases)
  }
}

# cumulative sum of counts associated with significant zones
sig_zones_cumsum = function(zones, y) {
  if (!is.list(zones)) stop("zones must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  vapply(zones, function(x) cumsum(y[x]), FUN.VALUE = numeric(1),
         USE.NAMES = FALSE)
}
