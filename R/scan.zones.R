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
  # number of regions
  N = nrow(coords)

  # compute intercentroid distance
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  
  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  mynn = nnpop(d, pop, ubpop)
  
  # use nearest neighbors to construct all zones
  zones = unlist(lapply(mynn, function(x) lapply(seq_along(x), function(i) x[seq_len(i)])), recursive = FALSE)

  # return only unique zones
  return(zones[distinct(zones)])
}

# argument checking for all scan tests
arg_check_scan_zones = 
  function(coords, pop, ubpop, longlat)
  {
    if(!(is.matrix(coords) | is.data.frame(coords))) stop("coords should be a matrix or a data frame")
    if(ncol(coords) != 2) stop("coords must have two columns")
    N = nrow(coords)
    if(length(pop) != N) stop("length(pop) != nrow(coords)")
    if(!is.numeric(pop)) stop("pop should be a numeric vector")
    if(length(ubpop) != 1 || !is.numeric(ubpop)) stop("ubpop should be a numeric vector of length 1")
    if(ubpop<= 0 || ubpop > 1) stop("ubpop should be a value between 0 and 1")
    if(length(longlat) != 1) stop("length(longlat) != 1")
    if(!is.logical(longlat)) stop("longlat should be a logical value")
  }


