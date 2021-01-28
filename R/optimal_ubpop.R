#' Optimal Population Upper Bound Statistics
#'
#' \code{optimal_ubpop} computes statistics for choosing an
#' optimal population upper bound. \code{ubpop_seq} is a sequence of values
#' to consider as the optimal choice of upper bound. The smallest value must
#' be at least \code{min(pop)/sum(pop)} and should generally be less than 0.5.
#'
#' @inheritParams scan.test
#' @param ubpop_seq A strictly increasing numeric vector with values between min(pop)/sum(pop) and 1. The default is \code{seq(0.01, 0.5, len = 50)}.
#'
#' @return Returns a \code{smerc_optimal_ubpop} object. This includes:
#' \item{ubpop_seq}{The sequence of population bounds considered}
#' \item{neg_lrt}{The negative likelihood ratio test statistics of \code{ubpop_seq} for the elbow method}
#' \item{gini_coef}{The Gini coefficients of \code{ubpop_seq} for the elbow method}
#' \item{eb_point}{A list with the index of the elbow point, the x-value of the elbow point, and the y-value of the elbow point}
#' \item{elbow_ubpop}{The population upperbound suggested by the elbow method}
#' \item{gini_ubpop}{The population upperbound suggested by the Gini method}
#' @seealso \code{\link{scan.test}}
#' @author Joshua French
#' @export
#' @references Blah
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' ubpop_stats = optimal_ubpop(coords = coords, cases  = nydf$cases,
#'                             pop = nydf$pop, nsim = 49,
#'                             ubpop = seq(0.05, 0.5, by = 0.05))
#' plot(ubpop_stats)
optimal_ubpop = function(coords, cases, pop,
                     ex = sum(cases) / sum(pop) * pop,
                     nsim = 499, alpha = 0.05,
                     ubpop_seq = seq(0.01, 0.5, len = 50),
                     longlat = FALSE, cl = NULL,
                     type = "poisson",
                     min.cases = 0,
                     simdist = "multinomial") {
  # argument checking
  type = match.arg(type, c("poisson", "binomial"))
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_optimal_ubpop(coords = coords, cases = cases,
                          pop = pop, ex = ex, nsim = nsim,
                          alpha = alpha, ubpop_seq = ubpop_seq,
                          longlat = longlat,
                          k = 1, w = diag(nrow(coords)),
                          type = type, simdist = simdist,
                          min.cases = min.cases)

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
  # logical vector of duplicate zones
  ldup = nndup(nn, N)

  # compute number of cases in each candidate zone
  # only keep non-duplicated zones
  yin = nn.cumsum(nn, cases)[!ldup]
  # compute number of cases in each candidate zone
  # only keep non-duplicated zones
  popin = nn.cumsum(nn, pop)[!ldup]
  # determine non-duplicated zones
  zones = nn2zones(nn)[!ldup]

  # compute test statistics for observed data
  if (type == "poisson") {
    ein = nn.cumsum(nn, ex)[!ldup]
    eout = sum(ex) - ein
    popout = NULL
    tobs = stat.poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    popout = tpop - popin
    tobs = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
  }
  # set tobs to 0 when min.cases restriction not satisfied
  tobs[yin < min.cases] = 0

  # determine tobs for each population upper bound
  # tobs_seq = tobs_keep_seq(tobs = tobs,
  #                          ubpop_seq = ubpop_seq,
  #                          popin = popin,
  #                          tpop = tpop,
  #                          yin = yin,
  #                          min.cases = min.cases)

  # get sequence of tobs filtered by population constraints
  tobs_seq = seq_filter_by_popin(x = tobs,
                                 popin = popin,
                                 tpop = tpop,
                                 ubpop_seq = ubpop_seq)

  # determine zones for each population upper bound
  # zones_seq = tobs_keep_seq(tobs = zones,
  #                           ubpop_seq = ubpop_seq,
  #                           popin = popin,
  #                           tpop = tpop,
  #                           yin = yin,
  #                           min.cases = min.cases)

  # get sequence of candidate filtered by population constraints
  zones_seq = seq_filter_by_popin(x = zones,
                                  popin = popin,
                                  tpop = tpop,
                                  ubpop_seq = ubpop_seq)

  # compute test statistics for simulated data for all candidate zones
  message("computing statistics for simulated data:")
  tall = seq_scan_sim(nsim = nsim, nn = nn, ty = ty,
                  ex = ex, type = type, ein = ein,
                  eout = eout, popin = popin,
                  popout = popout, tpop = tpop, cl = cl,
                  simdist = simdist, pop = pop,
                  min.cases = min.cases,
                  ldup = ldup,
                  return_type = "all")

  # for each population upper bound, determine the maximum test
  # statistic for each simulated data set
  message("partitioning statistics by ubpop_seq")
  # tsim_seq = tsim_ubpop_seq(ubpop_seq = ubpop_seq,
  #                           tall_yin_sim = tall_yin_sim,
  #                           min.cases = min.cases,
  #                           popin = popin,
  #                           tpop = tpop)

  tsim_seq = tsim_ubpop_seq2(ubpop_seq = ubpop_seq, tall = tall, popin = popin, tpop = tpop)

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

  # compute elbow point
  eb_point = elbow_point(ubpop_seq, neg_lrt)

  # for the significant clusters, compute the total cases and total ex
  # get sig pvalues for each population upper bound
  sig_pvalue_seq = lget(pruned_seq, name = "pvalue")
  # get sig zones for each population upper bound
  sig_zones_seq = lget(pruned_seq, name = "zones")

  # get number of cases in each sig zone for each population upper bound
  sig_yin_seq = lapply(sig_zones_seq, zones.sum, y = cases)
  # get expected cases in each sig zone for each population upper bound
  sig_ein_seq = lapply(sig_zones_seq, zones.sum, y = ex)

  # compute gini coefficients for sequence of population upper bounds
  gini_seq = mapply(gini_coeff,
                    casein = sig_yin_seq,
                    exin = sig_ein_seq,
                    MoreArgs = list(ty = ty))

  # assign any non-significant sums (in case only a MLC returned)
  all_sig = (sapply(pvalue_seq, min) < alpha)

  # compute gini coefficient
  gini_coef = gini_seq * all_sig

  structure(list(ubpop_seq = ubpop_seq,
                 neg_lrt = neg_lrt,
                 gini_coef = gini_seq * all_sig,
                 eb_point = eb_point,
                 elbow_ubpop = eb_point$x,
                 gini_ubpop = ubpop_seq[which.max(gini_coef)]),
            class = "smerc_optimal_ubpop")
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

# compute gini coefficient based on cases in each cluster
# expected cases in each cluster
# total number of cases
gini_coeff = function(casein, exin, ty) {
  # order by ex
  o = order(exin, decreasing = FALSE)

  # compute cumulative proportion of cases, ex in sig clusters
  # in order of ex
  cp_casein = cumsum(casein[o]) / ty
  cp_exin = cumsum(exin[o]) / ty

  # compute gini coefficient
  2 * (.5 - MESS::auc(c(0, cp_casein, 1), c(0, cp_exin, 1)))
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

# extract the values of x that satisfy the constraint
# that prop_popin (= popin/tpop) <= ubpop
filter_by_popin = function(x, prop_popin, ubpop) {
  # return the observed statistics that satisfy population constraint
  x[prop_popin <= ubpop]
}

# return a list with the sequence of x's that satisfy the constraint
# that prop_popin (= popin/tpop) <= ubpop
# for each ubpop in ubpop_seq
seq_filter_by_popin = function(x, popin, tpop, ubpop_seq) {
  if (length(x) != length(popin)) {
    stop("popin must have the same length as x")
  }
  prop_popin = popin/tpop
  lapply(seq_along(ubpop_seq), function(i) {
    filter_by_popin(x = x,
                    prop_popin = prop_popin,
                    ubpop = ubpop_seq[i])
  })
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


# based on a specific ubpop, determine the maximum test statistic
# for each simulated data set while accounting for the min.cases
# requirement
tsim_ubpop2 = function(ubpop, tall, prop_popin) {
  # determine which candidate zones satisfy population constraint
  keep_pop = (prop_popin <= ubpop)
  # for each simulated data set, determine the maximum test
  # statistic for the zones satisfying the population constraint
  # and the min.cases constraint
  sapply(seq_along(tall), function(i) {
    # max of statistic satisfying both constraints
    max(tall[[i]][keep_pop])
  })
}

# based on the sequence ubpop_seq, determine the maximum test statistic
# for each simulated data set while accounting for the min.cases
# requirement
tsim_ubpop_seq2 = function(ubpop_seq, tall, popin, tpop) {
  # for each ubpop in ubpop_seq, compute the maximum test statistic
  # for each simulated data set
  prop_popin = popin/tpop
  pbapply::pbsapply(X = ubpop_seq,
                    FUN = tsim_ubpop2,
                    tall = tall,
                    prop_popin = prop_popin, simplify = FALSE)
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
#' @param ubpop_seq An strictly increasing sequence of values between min(pop)/sum(pop) and 1.
#' @param longlat A logical. TRUE is great circle distance.
#' @param parallel Not used.
#' @param k Number of nearest neighbors. Not always needed.
#' @param w A spatial proximity matrix
#' @param type Statistic type
#' @param simdist Distribution of simulation
#' @param min.cases Minimum number of cases. Only for scan.test.
#' @return NULL
#' @noRd
arg_check_optimal_ubpop = function(coords, cases, pop, ex, nsim, alpha,
                                   ubpop_seq, longlat, k, w, type, simdist,
                                   min.cases) {
    arg_check_scan_test(coords = coords, cases = cases,
                        pop = pop, ex = ex, nsim = nsim,
                        alpha = alpha, ubpop = 0.1,
                        longlat = longlat,
                        k = 1, w = diag(nrow(coords)),
                        type = type, simdist = simdist,
                        min.cases = min.cases)
  lb = min(pop)/sum(pop)
  arg_check_ubpop_seq(ubpop_seq, lb)
}

# cumulative sum of counts associated with significant zones
sig_zones_cumsum = function(zones, y) {
  if (!is.list(zones)) stop("zones must be a list")
  if (!is.numeric(y)) stop("y must be a numeric vector")
  vapply(zones, function(x) cumsum(y[x]), FUN.VALUE = numeric(1),
         USE.NAMES = FALSE)
}
