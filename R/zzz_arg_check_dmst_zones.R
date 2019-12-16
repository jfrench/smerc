#' Argument checking for dmst.zones, dc.zones, and
#' mlink.zones functions
#'
#' Check the arguments of the dmst.zones, dc.zones, and
#' mlink.zones functions.
#' @return NULL
#' @export
#' @keywords internal
arg_check_dmst_zones = function(coords, cases, pop, w, ex,
                                ubpop, ubd, longlat, type,
                                progress = FALSE) {
  if (ncol(coords) != 2) {
    stop("coords should have 2 columns")
  }
  if (nrow(coords) != length(cases)) {
    stop("nrow(coords) != length(cases)")
  }
  if (length(cases) != length(pop)) {
    stop("length(cases) != length(pop)")
  }
  if (!is.numeric(cases)) {
    stop("cases must be numeric")
  }
  if (!is.numeric(pop)) {
    stop("pop must be numeric")
  }
  if (length(cases) != nrow(w)) {
    stop("length(cases) != nrow(w)")
  }
  if (length(cases) != length(ex)) {
    stop("length(cases) != length(ex)")
  }
  if (length(ubpop) != 1 | !is.numeric(ubpop)) {
    stop("ubpop should be a single number")
  }
  if (ubpop <= 0 | ubpop > 1) {
    stop("ubpop not in (0, 1]")
  }
  if (length(ubd) != 1 | !is.numeric(ubd)) {
    stop("ubd should be a single number")
  }
  if (ubd <= 0 | ubd > 1) {
    stop("ubd not in (0, 1]")
  }
  if (length(longlat) != 1 || !is.logical(longlat)) {
    stop("longlat must be a single logical value")
  }
  if (!is.element(type, c("maxonly", "pruned", "all"))) {
    stop("Invalid type")
  }
  if (length(progress) != 1 | !is.logical(progress)) {
    stop("progress must be a logical value")
  }
}
