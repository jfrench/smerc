#' Compute Poisson test statistic
#'
#' Efficiently compute a vector of Poisson test statistics.
#' This function does no sanity checking. We encourage you
#' to compare the results against
#' \code{\link[smerc]{scan.stat}} for correctness.
#' \code{stat_poisson_adj} is a C++ version implementation
#' of the code and should be faster.
#' \code{stat_binomial_adj} is a C++ version implementation
#' of \code{\link[smerc]{stat.binom}}.
#'
#' @inheritParams scan.stat
#' @inheritParams scan.test
#' @param logein The \code{log} of the expected number of
#'   cases in each candidate zone.
#' @param logeout The \code{log} of the expected number of
#'   cases outside of each candidate zone.
#' @param pen The eccentricity penalty for each candidate
#'   zone.
#' @param logpopin The \code{log} of the population in each
#' candidate zone.
#' @param logpopout The \code{log} of the population outside
#' of each candidate zone.
#' @param return.max A logical value. Default is \code{FALSE}. If
#' \code{TRUE}, then only the maximum statistic is returned.
#'
#' @return A numeric vector.
#' @keywords internal
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' enn = elliptic.nn(coords, nydf$pop, ubpop = 0.015)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' yin = nn.cumsum(enn$nn, cases)
#' ein = nn.cumsum(enn$nn, ex)
#' logein = log(ein)
#' logeout = log(ty - ein)
#' pen = elliptic.penalty(a = 0.5, enn$shape_all)
#' stat.poisson.adj(yin, ty, logein, logeout,
#'                  a = 0.5, pen = pen, return.max = TRUE)
stat.poisson.adj = function(yin, ty, logein, logeout, a = 0,
                            pen = 1, min.cases = 2,
                            return.max = FALSE) {
  yout = ty - yin
  # determine if there will be any problematic statistics
  good = which(yin >= min.cases)
  # create vector for storage
  tall = numeric(length(yin))
  # log ratio observed/expected (in and out) for good locations
  lrin =  log(yin[good]) - logein[good]
  lrout = log(yout[good]) - logeout[good]
  # compute statistics for good locations
  tall[good] = yin[good] * lrin + yout[good] * lrout
  # if indicator not satisfied (yin/ein > yout/eout), set to 0
  tall[good][lrin < lrout] = 0
  if (a > 0) {
    tall = tall * pen
  }
  if (return.max) {
    return(max(tall))
  } else {
    return(tall)
  }
}

#' @export
#' @rdname stat.poisson.adj
stat_poisson_adj = function(yin, ty, logein, logeout, a = 0,
                            pen = 1, min.cases = 2,
                            return.max = FALSE) {
  .Call(`_smerc_stat_poisson_adj_cpp`, yin, ty, logein, logeout, a, pen, min.cases, return.max)
}

#' @export
#' @rdname stat.poisson.adj
stat_binom_adj <- function(yin, ty, popin, popout,
                           logpopin, logpopout, tpop,
                           min.cases = 2,
                           return.max = FALSE) {
  .Call(`_smerc_stat_binom_adj_cpp`, yin, ty, popin, popout,
        logpopin, logpopout, tpop, min.cases, return.max)
}
