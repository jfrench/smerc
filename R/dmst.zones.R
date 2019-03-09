#' Determine zones for the Dynamic Minimum Spanning Tree scan test
#'
#' \code{dmst.zones} determines the zones for the Dynamic
#' Minimum Spanning Tree scan test (\code{\link{dmst.test}}).  The
#' function returns the zones, as well as the associated
#' test statistic, cases in each zone, the expected number
#' of cases in each zone, and the population in each zone.
#'
#' Every zone considered must have a total population less
#' than \code{ubpop * sum(pop)}.  Additionally, the maximum
#' intercentroid distance for the regions within a zone must
#' be no more than \code{ubd * the maximum intercentroid
#' distance across all regions}.
#' @inheritParams dmst.test
#' @inheritParams mst.all
#' @inheritParams flex.zones
#' @inherit dc.zones return
#' @author Joshua French
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and
#'   Neto, S.J.F. (2006). Fast detection of arbitrarily
#'   shaped disease clusters, Statistics in Medicine, 25,
#'   723-742.  <https://doi.org/10.1002/sim.2411>
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' # find zone with max statistic starting from each individual region
#' all_zones = dmst.zones(coords, cases = floor(nydf$cases),
#'                         nydf$pop, w = nyw, ubpop = 0.25,
#'                         ubd = .25, longlat = TRUE)
dmst.zones = function(coords, cases, pop, w,
                       ex = sum(cases) / sum(pop) * pop,
                       ubpop = 0.5, ubd = 1, longlat = FALSE,
                       cl = NULL, progress = TRUE) {
  # sanity checking
  arg_check_dmst_zones(coords = coords, cases = cases,
                       pop = pop, w = w, ex = ex,
                       ubpop = ubpop, ubd = ubd,
                       longlat = longlat, type = "all",
                       progress = progress)
  # setup various arguments and such
  ty = sum(cases)   # total number of cases
  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # find all neighbors from each starting zone within distance upperbound
  nn = nndist(d, ubd)
  # get zones and relevant information
  out = mst.all(neighbors = nn, cases = cases, pop = pop,
                w = w, ex = ex, ty = ty, max_pop = max_pop,
                type = "all", nlinks = "one", early = FALSE,
                cl = cl, progress = progress)
  prep.mst(out)
}

#' Return nicely formatted results from mst.all
#'
#' Return nicely formatted results from mst.all
#' @return NULL
#' @export
#' @keywords internal
prep.mst = function(mstout) {
  nn = lapply(mstout, getElement, name = "locids")
  loglikrat = unlist(lapply(mstout, getElement,
                            name = "loglikrat"))
  cases = unlist(lapply(mstout, getElement, name = "cases"))
  expected = unlist(lapply(mstout, getElement,
                           name = "expected"))
  population = unlist(lapply(mstout, getElement,
                             name = "population"))
  return(list(zones = nn2zones(nn),
              loglikrat = loglikrat,
              cases = cases,
              expected = expected,
              population = population))
}
