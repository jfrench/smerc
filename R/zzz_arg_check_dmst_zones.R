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
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_w(w, N)
  arg_check_ex(ex, N)
  arg_check_ubpop(ubpop)
  arg_check_ubd(ubd)
  arg_check_longlat(longlat)
  if (!is.element(type, c("maxonly", "pruned", "all"))) {
    stop("Invalid type")
  }
  if (length(progress) != 1 | !is.logical(progress)) {
    stop("progress must be a logical value")
  }
}
