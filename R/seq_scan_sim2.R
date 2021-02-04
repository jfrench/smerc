#' Perform \code{scan.test} on simulated data sequentially
#'
#' \code{seq_scan_sim2} efficiently performs \code{\link{scan.test}} on a
#' simulated data set.  The function is meant to be used internally by the
#' \code{\link{optimal_ubpop}} function.
#'
#' @inheritParams scan.sim
#' @inheritParams scan.test
#' @param ldup A logical vector indicating positions of duplicated zones. Not
#'   intended for user use.
#' @param lseq_zones A list of logical vectors specifying the sequence of
#'   relevant zones based on ubpop constraints
#' @return If a matrix with the maximum test statistics. Each row of the matrix
#'   provides the maximum test statistics for a particular set of population
#'   constraints.
#' @export
#' @keywords internal
seq_scan_sim2 = function(nsim = 1, nn, ty, ex, type = "poisson",
                    ein = NULL, eout = NULL,
                    tpop = NULL, popin = NULL, popout = NULL,
                    cl = NULL,
                    simdist = "multinomial",
                    pop = NULL,
                    min.cases = 0,
                    ldup = NULL,
                    lseq_zones) {
  # match simdist with options
  simdist = match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_seq_scan_sim2(nsim = nsim, ty = ty, ex = ex, type = type,
                nn = nn, ein = ein, eout = eout, tpop = tpop,
                popin = popin, popout = popout, static = TRUE,
                simdist = simdist, pop = pop,
                w = diag(length(ex)),
                ldup = ldup,
                lseq_zones = lseq_zones)
  # assume there are no duplicates if ldup not provided
  if(is.null(ldup)) {
    ldup = rep(FALSE, length(unlist(nn)))
  }

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pbsapply(seq_len(nsim), function(i) {
    # simulate new data
    if (simdist == "multinomial") {
      ysim = stats::rmultinom(1, size = ty, prob = ex)
    } else if (simdist == "poisson") {
      ysim = stats::rpois(length(ex), lambda = ex)
      ty = sum(ysim)
      mult = ty / sum(ex)
      ein = ein * mult
      eout = eout * mult
    } else if (simdist == "binomial") {
      ysim = stats::rbinom(n = length(ex), size = pop,
                           prob = ex / pop)
      ty = sum(ysim)
    }
    # compute test statistics for each zone that aren't duplicated
    yin = nn.cumsum(nn, ysim)[!ldup]
    if (type == "poisson") {
      tall = stat.poisson(yin, ty - yin, ein, eout)
    } else if (type == "binomial") {
      tall = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
    }
    tall[yin < min.cases] = 0

    # return sequence of maximum statistic for each element
    # of lseq_zones
    sapply(lseq_zones, function(lzones) {
      max(tall[lzones])
    })
  })
  return(tsim)
}

#' Argument checking for *.sim functions
#'
#' Check the arguments of the \code{*.sim} functions.
#'
#' @param nsim Number of simulations
#' @param ty Total number of cases
#' @param ex Expected counts
#' @param type Type of statistic
#' @param nn List of nn (e.g., nnpop function)
#' @param zones List of zones (e.g., scan.zones)
#' @param ein List of expected in each zone
#' @param eout List of expected out of each zone
#' @param tpop Total population
#' @param popin Population in each zone
#' @param popout Population outside of each zone
#' @param w Spatial adjacency matrix
#' @param pop Vector of populations
#' @param ubpop Population upperbound
#' @param static Static zones. Logical. TRUE for scan.test.
#' FALSE for uls.test.
#' @param simdist Simulation distribution.
#' @param ldup A logical vector indicating positions of duplicated zones. Not intended for user use.
#' @param lseq_zones A list of logical vectors specifying the sequence of relevant zones based on ubpop constraints

#' @return NULL
#' @noRd
arg_check_seq_scan_sim2 = function(nsim, ty, ex, type,
                         nn = NULL, zones = NULL,
                         ein = NULL, eout = NULL,
                         tpop = NULL, popin = NULL,
                         popout = NULL, w = NULL,
                         pop = NULL, ubpop = NULL,
                         static = FALSE,
                         simdist = "multinomial",
                         ldup = NULL,
                         lseq_zones) {
  arg_check_nsim(nsim)
  arg_check_ty(ty)
  N = length(ex)
  arg_check_ex(ex, N)
  arg_check_type(type)
  if (!is.null(nn)) {
    if (!is.list(nn)) stop("nn must be a list")
    nz = sum(sapply(nn, length))
  }
  if (!is.null(zones)) {
    if (!is.list(zones)) stop("zones must be a list")
  }
  arg_check_w(w, N)
  if (!is.null(ubpop)) {
    arg_check_ubpop
  }
  if (!is.null(pop)) {
    arg_check_pop(pop, N)
  }
  arg_check_simdist(simdist)
  if (simdist == "binomial" & is.null(pop)) {
    stop("pop must be specified when simdist == 'binomial'")
  }
  arg_check_ldup(ldup)
  if (length(ldup) != length(unlist(nn))) {
    stop("length(ldup) doesn't match the number of total candidate zones (length(unlist(nn)))")
  }
  if (!is.list(lseq_zones)) {
    stop("lseq_zones must be a list")
  }
  lseq_zones_lengths = sapply(lseq_zones, length)
  if (min(lseq_zones_lengths) != max(lseq_zones_lengths)) {
    stop("lseq_zones_lengths must all be the same")
  }
  # each element of lseq_zones should have length equal to the
  # number of non-duplicated elements (sum(!ldup))
  if (min(lseq_zones_lengths) != sum(!ldup)) {
    stop("sum(!ldup) should be the same as the length of each element of lseq_zones")
  }
}

