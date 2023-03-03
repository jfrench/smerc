#' Compute elliptic penalty
#'
#' Compute eccentricity penalty for elliptic scan method.
#'
#' @param a Penalty scale
#' @param shape Shape of ellipse.
#'
#' @return A vector of penalities
#' @export
#'
#' @examples
#' elliptic.penalty(a = 0.5, shape = c(1, 1.5, 2))
elliptic.penalty <- function(a, shape) {
  if (length(a) != 1 | !is.numeric(a) | min(a) < 0) {
    stop("a must be a single positive numeric value")
  }
  if (!is.numeric(a) | min(shape) < 1) {
    stop("shape values must be at least 1")
  }
  (4 * shape / (shape + 1)^2)^a
}
