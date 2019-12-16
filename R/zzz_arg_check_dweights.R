#' Argument checking for dweights/tango.weights
#'
#' @param coords Matrix of coordinates
#' @param kappa Scale parameter for weight function
#' @param longlat Logical value indicating whether distance
#' should be Euclidean (FALSE) or great citcle (TRUE)
#' @param type Type of weights (basic, rogerson, or tango)
#' @param pop vector of population values
#'
#' @return NULL
#' @export
#'
#' @keywords internal
arg_check_dweights = function(coords, kappa, longlat, type, pop) {
  if (!(is.matrix(coords) | is.data.frame(coords))) {
    stop("coords should be a matrix or a data frame")
  }
  if (ncol(coords) != 2) {
    stop("coords must have two columns")
  }
  N = nrow(coords)
  if (length(kappa) != 1 || !is.numeric(kappa)) {
    stop("kappa should be a numeric vector of length 1")
  }
  if (kappa <= 0) stop("kappa must be positive")
  if (length(longlat) != 1) stop("length(longlat) != 1")
  if (!is.logical(longlat)) {
    stop("longlat should be a logical value")
  }
  if (length(type) != 1) {
    stop("type must be a single character")
  }
  if (!is.element(type, c("basic", "rogerson", "tango"))) {
    stop("invalid type")
  }
  if (!is.null(pop)) {
    if (length(pop) != N) {
      stop("length(pop) != nrow(coords)")
    }
    if (!is.numeric(pop)) {
      stop("pop should be a numeric vector")
    }
  }
  if (type == "rogerson") {
    if (is.null(pop)) {
      stop("pop must be provided when type = 'rogerson'")
    }
    if (any(pop == 0)) {
      stop("regions cannot contain zero population when type = 'rogerson'")
    }
  }
}
