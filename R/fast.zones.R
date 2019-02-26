#' Determine sequence of fast subset scan zones
#'
#' \code{fast.zones} determines the unique zones obtained by
#' implementing the fast subset scan method of Neill (2012).
#'
#' The zones returned must have a total population less than
#' \code{ubpop * sum(pop)} of all regions in the study area.
#'
#' @inheritParams uls.test
#' @return Returns a vector of regions to sequentially and
#'   cumulatively consider for clustering.
#' @author Joshua French
#' @export
#' @references Neill, D. B. (2012), Fast subset scan for
#'   spatial pattern detection. Journal of the Royal
#'   Statistical Society: Series B (Statistical
#'   Methodology), 74: 337-360.
#'   <doi:10.1111/j.1467-9868.2011.01014.x>
#' @examples
#' data(nydf)
#' fast.zones(cases = nydf$cases, pop = nydf$population)
fast.zones = function(cases, pop, ubpop = 0.5) {
  if (length(cases) != length(pop)) stop('length(cases) != length(pop)')
  if (length(ubpop) != 1 | !is.numeric(ubpop)) stop("ubpop should be a single number")
  if (ubpop <= 0 | ubpop > 1) stop("ubpop not in (0, 1]")
  
  # order rates from largest to smallest
  or = order(cases/pop, decreasing = TRUE)
  #
  max_pop = sum(pop) * ubpop
  return(or[cumsum(pop[or]) <= max_pop])
}