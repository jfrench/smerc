#' Perform \code{scan.test} on simulated data
#'
#' \code{scan.sim} efficiently performs
#' \code{\link{scan.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{scan.test}} function, but is informative for
#' better understanding the implementation of the test.
#' 
#' @inheritParams flex.sim
#' @inheritParams nn.cumsum
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' nn = nnpop(d, pop = nydf$pop, ubpop = 0.1)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' yin = nn.cumsum(nn, cases)
#' ein = nn.cumsum(nn, ex)
#' tsim = scan.sim(nsim = 2, nn, ty, ex, ein = ein, eout = ty - ein)
scan.sim = function(nsim = 1, nn, ty, ex, type = "poisson", 
                    ein = NULL, eout = NULL, 
                    tpop = NULL, popin = NULL, popout = NULL, 
                    cl = NULL) {
  
  arg_check_scan_sim(nsim, nn, ty, ex, ein, eout,
                     tpop, popin, popout, type)
  
  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # compute test statistics for each zone
    yin = nn.cumsum(nn, ysim)
    if (type == "poisson") {
      tall = stat.poisson(yin, ty - yin, ein, eout)
    } else if (type == "binomial") {
      tall = stat.binom(yin, ty - yin, ty, 
                        popin, popout, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}

arg_check_scan_sim = function(nsim, nn, ty, ex, ein, eout,
                              tpop, popin, popout, type) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  if (!is.list(nn)) stop("nn must be a list")
  if (length(ty) != 1) stop("ty must be a single number")
  if (!is.numeric(ex)) stop("ex must be a numeric vector")
  if (length(type) != 1 | !is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
  nz = sum(sapply(nn, length))
  if (type == "poisson") {
    if (is.null(ein) | is.null(eout)) {
      stop("ein and eout must be provided when type is 'poisson'")
    }
    if (nz != length(ein)) {
      stop("length(ein) has wrong length")
    }
    if (length(ein) != length(eout)) {
      stop("length(ein) != length(eout)")
    }
  }
  if (type == "binomial") {
    if (is.null(popin) | is.null(popout) | is.null(tpop)) {
      stop("popin, popout, and tpop must be provided when type is 'binomial'")
    }
    if (nz != length(popin)) {
      stop("length(popin) has wrong length")
    }
    if (nz != length(popout)) {
      stop("length(popin) != length(popout)")
    }
    if (length(tpop) != 1) {
      stop("tpop must be a single number")
    }
  }
}
