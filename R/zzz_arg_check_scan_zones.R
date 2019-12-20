#' Argument checking for scan.zones
#'
#' Check the arguments of the scan.zones function
#' @param coords A matrix of coordinates
#' @param pop A vector of population for each centroid
#' @param ubpop The population upperbound
#' @param longlat TRUE means great circle distance
#' @return NULL
#' @noRd
# argument checking for all scan zones
arg_check_scan_zones =
  function(coords, pop, ubpop, longlat) {
    arg_check_coords(coords)
    N = nrow(coords)
    arg_check_pop(pop, N)
    arg_check_ubpop(ubpop)
    arg_check_longlat(longlat)
  }
