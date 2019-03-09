#' Determine zones for the spatial scan test
#'
#' \code{scan.zones} determines the unique candidate
#' zones to consider for the circular spatial scan test of
#' Kulldorff (1997).
#'
#' @inheritParams scan.test
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial
#' scan statistic. Communications in Statistics - Theory and
#' Methods, 26(6): 1481-1496,
#' <doi:10.1080/03610929708831995>
#' @examples
#' data(nydf)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' zones = scan.zones(coords = coords, pop = nydf$pop,
#'                    ubpop = 0.1, longlat = TRUE)
scan.zones = function(coords, pop, ubpop = 0.5, longlat = FALSE) {
  # argument checking
  arg_check_scan_zones(coords, pop, ubpop, longlat)

  # compute intercentroid distance
  d = sp::spDists(as.matrix(coords), longlat = longlat)

  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = nnpop(d, pop, ubpop)

  # use nearest neighbors to construct all zones
  zones = nn2zones(nn)

  # return only unique zones
  return(zones[distinct(zones)])
}
