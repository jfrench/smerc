#' Perform \code{cepp.test} on simulated data
#'
#' \code{cepp.sim} efficiently performs
#' \code{\link{cepp.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{cepp.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @inheritParams flex.sim
#' @inheritParams cepp.test
#' @param nn A list of nearest neighbors produced by \code{\link{casewin}}.
#' @param w A list that has the weights associated with each
#' region of each element of \code{nn}.
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' d = sp::spDists(as.matrix(coords), longlat = TRUE)
#' nn = casewin(d, cases = nydf$pop, cstar = 15000)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' # find smallest windows with at least n* pop
#' nstar = 15000
#' cwins = casewin(d, cases = nydf$pop, cstar = nstar)
#' # number of regions in each window
#' l = sapply(cwins, length)
#'
#' # determine population in each window
#' pop_cwins = nn.cumsum(cwins, nydf$pop, simplify = FALSE)
#' # determine weights for each region in each window
#' w = lapply(seq_along(l), function(i) {
#'   out = rep(1, l[i])
#'   out[l[i]] = 1 - pop_cwins[[i]][l[i] - 1]/nstar
#'   return(out)
#' })
#' tsim = cepp.sim(1, nn = cwins, ty = ty, ex = ex, w = w)
cepp.sim = function(nsim = 1, nn, ty, ex, w, simtype = "multinomial") {

  # number of regions in each window
  l = sapply(nn, length)

  # compute max test stat for nsim simulated data sets
  tsim = pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    if (simtype == "multinomial") {
      ysim = stats::rmultinom(1, size = ty, prob = ex)
    } else {
      ysim = stats::rpois(length(ex), lambda = ex)
    }
    cstar = sapply(seq_along(l), function(i) {
      sum(ysim[nn[[i]]] * w[[i]])
    })
    max(cstar)
  })
  unlist(tsim, use.names = FALSE)
}
