#' Argument checking for dweights/tango.weights
#'
#' @param coords Matrix of coordinates
#' @param kappa Scale parameter for weight function
#' @param longlat Logical value indicating whether distance
#' should be Euclidean (FALSE) or great circle (TRUE)
#' @param type Type of weights (basic, rogerson, or tango)
#' @param pop vector of population values
#'
#' @return NULL
#' @noRd
arg_check_dweights = function(coords, kappa, longlat, type, pop) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_dweights_kappa(kappa)
  arg_check_longlat(longlat)
  if (!is.null(pop)) {
    arg_check_pop(pop, N)
  }
  arg_check_dweights_type(type)
  if (type == "rogerson" & is.null(pop)) {
    stop("pop must be provided when type = 'rogerson'")
  }
}
