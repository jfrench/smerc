#' Perform \code{scan.test} on simualated data
#'
#' \code{scan.sim} efficiently performs
#' \code{\link{scan.test}} on a simulated data set.  The
#' function is meant to be used internally by the
#' \code{\link{scan.test}} function, but is informative for
#' better understanding the implementation of the test.
#' 
#' \code{scan.sim} and \code{flex.sim} are identicial, as
#' the only thing that changes is the \code{zones} considered.
#'
#' @inheritParams flex.sim
#'
#' @return A vector with the maximum test statistic for each
#'   simulated data set.
#' @export
#'
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' zones = scan.zones(coords, pop = nydf$pop, longlat = TRUE, ubpop = 0.1)
#' cases = floor(nydf$cases)
#' ty = sum(cases)
#' ex = ty/sum(nydf$pop) * nydf$pop
#' yin = sapply(zones, function(x) sum(cases[x]))
#' ein = sapply(zones, function(x) sum(ex[x]))
#' tsim = scan.sim(nsim = 2, zones, ty, ex, ein = ein, eout = ty - ein)
scan.sim = flex.sim