#' Perform \code{scan.test} on simulated data
#'
#' \code{scan.sim} efficiently performs
#' \code{\link{scan.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{scan.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams flex.sim
#' @param nn A list of nearest neighbors produced by \code{\link{nnpop}}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' nn = scan.nn(d, pop = nydf$pop, ubpop = 0.1)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' yin = nn.cumsum(nn, cases)
#' ein = nn.cumsum(nn, ex)
#' tsim = scan.sim(nsim = 1, nn, ty, ex, ein = ein, eout = ty - ein)
scan.sim = function(nsim = 1, nn, ty, ex, type = "poisson",
                    ein = NULL, eout = NULL,
                    tpop = NULL, popin = NULL, popout = NULL,
                    cl = NULL) {
  arg_check_sim(nsim = nsim, ty = ty, ex = ex, type = type,
                nn = nn, ein = ein, eout = eout, tpop = tpop,
                popin = popin, popout = popout, static = TRUE)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # compute test statistics for each zone
    yin = nn.cumsum(nn, ysim)
    if (type == "poisson") {
      tall = stat.poisson(yin, ty - yin, ein, eout)
    } else if (type == "binomial") {
      tall = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
    }
    max(tall)
  })
  unlist(tsim, use.names = FALSE)
}
