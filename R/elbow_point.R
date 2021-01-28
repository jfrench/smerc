#' Compute Elbow Point
#'
#' \code{elbow_point} computes the elbow point based on the maximum distance
#' between each point and the line passing through the end points.
#'
#' @param x A numeric vector
#' @param y A numeric vector
#' @author Joshua French and Mohammad Meysami
#' @return A list with the index (\code{idx}), x-value (\code{x}) and
#' y-value (\code{y}) of the elbow point.
#' @export
#' @references  \url{https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line}
#' @examples
#' # generate some data
#' x = c(0, 0.5, 1)
#' y = c(1, 0.1, 0)
#' # plot data (the second point is clearly the elbow)
#' plot(x, y)
#' elbow_point(x, y)
elbow_point = function(x, y) {
  arg_check_elbow_point(x = x, y = y)

  # start and end point positions
  wmin_x = which.min(x)
  wmax_x = which.max(x)

  # start point
  x_start = x[wmin_x]
  y_start = y[wmin_x]
  # end point
  x_end = x[wmax_x]
  y_end = y[wmax_x]

  # distance between x, y start and end
  dx = x_end - x_start
  dy = y_end - y_start
  d_start_end = sqrt(dx^2 + dy^2)
  # distance between each point in (x, y) and the line passing through the
  # points (x_start, y_start) and (x_end, y_end)
  d = abs(dx * (y_start - y) - (x_start - x) * dy) / d_start_end

  wmax_d = which.max(d)
  return(list(idx = wmax_d, x = x[wmax_d], y = y[wmax_d]))
}

#' Check arguments of elbow_point
#'
#' @param x A numeric vector
#' @param y A numeric vector (same length as x)
#'
#' @return NULL
#' @noRd
arg_check_elbow_point = function(x, y) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  if (!is.numeric(x)) {
    stop("x must be a vector")
  }
  if (!is.numeric(y)) {
    stop("y must be numeric")
  }
  if (!is.numeric(y)) {
    stop("y must be a vector")
  }
  if (length(x) != length(y)) {
    stop("x and y must have the same length")
  }
  if (length(x) < 2) {
    stop("x and y must have at least 2 points")
  }
  if (length(x) < 3) {
    warning("x and y should have at least 3 points to find the elbow")
  }
  if (anyDuplicated(x)) {
    stop("")
  }
}
