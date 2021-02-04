#' Sequential Scan Test
#'
#' Performs a series of sequential scan tests by ubpop. Only minimal information is returned
#' for each test.
#'
#' @inheritParams optimal_ubpop
#'
#' @return Returns a list of \code{length(ubpop_seq)}. Each element of the
#' list has "pruned" set of results (as a list) for a scan test. Element of
#' the pruned list has the elements \code{tobs}, \code{zones}, \code{pvalue}.
#' @seealso \code{\link{scan.test}}
#' @author Joshua French
#' @export
#' @keywords internal
seq_scan_test = function(coords, cases, pop,
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
  arg_check_seq_scan_test(coords = coords, cases = cases,
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

  # compute logical vector specifying sequence of zones based
  # on ubpop constraints
  lseq_zones = lapply(ubpop_seq, function(ubpop) {
    popin/tpop <= ubpop}
  )

  # get sequence of tobs filtered by population constraints
  tobs_seq = lapply(lseq_zones, function(lseq) { tobs[lseq] })

  # get sequence of candidate zones filtered by population constraints
  zones_seq = lapply(lseq_zones, function(lseq) { zones[lseq] })

  # compute maximum test statistics for simulated data for all candidate zones
  # for each population constraint
  # each row contains the set of maximum statistics for each population constraint
  message("computing statistics for simulated data:")
  tsim_seq = seq_scan_sim2(nsim = nsim, nn = nn, ty = ty,
                           ex = ex, type = type, ein = ein,
                           eout = eout, popin = popin,
                           popout = popout, tpop = tpop, cl = cl,
                           simdist = simdist, pop = pop,
                           min.cases = min.cases,
                           ldup = ldup,
                           lseq_zones = lseq_zones)

  message("computing p-values")
  pvalue_seq = pbapply::pblapply(seq_along(ubpop_seq), function(i) {
    mc_pvalue_cpp(tobs = tobs_seq[[i]], tsim = tsim_seq[i, ])
  })

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
  return(pruned_seq)
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
arg_check_seq_scan_test = function(coords, cases, pop, ex, nsim, alpha,
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
