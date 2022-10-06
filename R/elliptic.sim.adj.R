#' Perform \code{elliptic.test} on simulated data
#'
#' \code{elliptic.sim} efficiently performs
#' \code{\link{elliptic.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{elliptic.test}} function, but is informative
#' for better understanding the implementation of the test.
#'
#' @inheritParams scan.sim
#' @param nn A list of nearest neighbors produced by
#'   \code{\link{elliptic.nn}}.
#' @inheritParams elliptic.test
#' @param shape_all A vector of the shapes associated with
#'   all of the possible zones constructed from \code{nn}.
#'   This can be obtained from \code{\link{elliptic.nn}}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' pop = nydf$pop
#' enn = elliptic.nn(coords, pop, ubpop = 0.5)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(pop) * pop
#' yin = nn.cumsum(enn$nn, cases)
#' ein = nn.cumsum(enn$nn, ex)
#' set.seed(1)
#' tsim = elliptic.sim(nsim = 3, nn = enn$nn, ty = ty, ex = ex,
#'                     a = 0.5, shape_all = enn$shape_all,
#'                     ein = ein, eout = ty - ein)
#' logein = log(ein)
#' logeout = log(ty - ein)
#' pen = elliptic.penalty(0.5, enn$shape_all)
#' set.seed(1)
#' tsim2 = elliptic.sim.adj(nsim = 3, ex = ex,
#'                          nn = enn$nn, ty = ty,
#'                          logein = logein, logeout = logeout,
#'                          a = 0.5, pen = pen)
#' all.equal(tsim, tsim2)
elliptic.sim.adj = function(nsim = 1, ex, nn, ty, logein, logeout,
                            a, pen, min.cases = 2, cl = NULL) {
  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # compute test statistics for each zone
    yin = nn.cumsum(nn, ysim)
    stat.poisson.adj(yin = yin,
                     ty = ty,
                     logein = logein, logeout = logeout,
                     a = a, pen = pen,
                     min.cases = min.cases,
                     return.max = TRUE)
  }, cl = cl)
  unlist(tsim, use.names = FALSE)
}

