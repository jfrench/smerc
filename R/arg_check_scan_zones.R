#' Argument checking for scan.zones
#'
#' Check the arguments of the scan.zones function
#' @return NULL
#' @export
#' @keywords internal
# argument checking for all scan zones
arg_check_scan_zones =
  function(coords, pop, ubpop, longlat) {
    if (!(is.matrix(coords) | is.data.frame(coords))) {
      stop("coords should be a matrix or a data frame")
    }
    if (ncol(coords) != 2) {
      stop("coords must have two columns")
    }
    N = nrow(coords)
    if (length(pop) != N) {
      stop("length(pop) != nrow(coords)")
    }
    if (!is.numeric(pop)) {
      stop("pop should be a numeric vector")
    }
    if (length(ubpop) != 1 || !is.numeric(ubpop)) {
      stop("ubpop should be a numeric vector of length 1")
    }
    if (ubpop <= 0 || ubpop > 1) {
      stop("ubpop should be a value between 0 and 1")
    }
    if (length(longlat) != 1) {
      stop("length(longlat) != 1")
    }
    if (!is.logical(longlat)) {
      stop("longlat should be a logical value")
    }
  }
