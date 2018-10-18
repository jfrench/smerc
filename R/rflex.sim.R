#' Perform \code{rflex.test} on simualated data
#'
#' \code{rflex.sim} efficiently performs
#' \code{\link{rflex.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{rflex.test}} function, but is informative for
#' better understanding the implementation of the test.
#'
#' @param nsim A positive integer indicating the number of
#'   simulations to perform.
#' @param nn A matrix of the k nearest neighbors for the
#'   regions described by \code{w}.
#' @inheritParams rflex.test
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' # determine knn
#' coords = with(nydf, cbind(longitude, latitude))
#' nn = knn(coords, longlat = TRUE, k = 50)
#' # determine expected number of cases in each region
#' cases = floor(nydf$cases)
#' pop = nydf$pop
#' ex = pop * sum(cases)/sum(pop)
#' tsim = rflex.sim(nsim = 5, nn = nn, w = nyw, ex = ex)
rflex.sim = function(nsim = 1, nn, w, ex, type = "poisson", 
                     alpha1 = 0.2, cl = NULL) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  ty = sum(ex)
  pbapply::pbsapply(seq_len(nsim), function(i) {
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    zones = rflex.zones(nn = nn, w = w, 
                        cases = ysim, ex = ex, 
                        alpha1 = alpha1,
                        cl = NULL, progress = FALSE)
    ein = unlist(lapply(zones, function(x) sum(ex[x])), use.names = FALSE)
    yin = unlist(lapply(zones, function(x) sum(ysim[x])))
    tall = scan.stat(yin, ein, ty - ein, ty, type)
    max(tall)
  })
}