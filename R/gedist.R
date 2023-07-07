#' Compute distance for geographic coordinates
#'
#' \code{gedist} computes the distance between the
#' coordinates in \code{x} and \code{y}. If \code{y} isn't
#' supplied, then the distances are computed between the
#' coordinates in \code{x} alone. Otherwise, the pairwise
#' distances between the points in \code{x} and \code{y} are
#' computed. If \code{longlat = TRUE}, then the great circle
#' distance is computed. \code{eucdist} is a simplified
#' version of \code{gedist} that computes Euclidean
#' distances alone while \code{gcdist} is a simplified
#' version of \code{gedist} that computes great circle
#' distance alone.
#'
#' The algorithm used when \code{longlat = TRUE} is a C++
#' port of the C code written by Roger Bivand for
#' \code{\link[sp]{spDists}}, which appears to be based on a
#' special case of the Vincenty formula with a slight
#' correction based on the WGS84 flattening constant. See
#' \url{https://en.wikipedia.org/wiki/Great-circle_distance}.
#'
#' @param x A two-dimensional matrix of coordinates.
#' @param y A two-dimensional matrix of coordinates.
#' @param longlat A logical value indicating whether
#'   Euclidean distance (\code{longlat = FALSE}) or great
#'   circle distance (\code{longlat = FALSE}) should be
#'   computed. The default is \code{longlat = FALSE}.
#'
#' @return A matrix of distances
#' @export
#' @examples
#' coords = matrix(runif(10), ncol = 2)
#' # euclidean distance
#' d = gedist(coords)
#' all.equal(d, as.matrix(dist(coords)),
#'           check.attributes = FALSE)
#' all.equal(gedist(coords), eucdist(coords))
#'
#' # great circle distance
#' all.equal(gedist(coords, longlat = TRUE),
#'           gcdist(coords))
gedist = function(x, y = NULL, longlat = FALSE) {
  arg_check_coords(x, "x")
  if (!is.null(y)) {
    arg_check_coords(y, "y")
  }

  # if only 1 set of coordinates
  if (is.null(y)) {
    if (!longlat) {
      # euclidean distance, 1 set of coordinates
      eucdist1(x = x[,1], y = x[,2], eps = .Machine$double.eps)
    } else {
      # great circle distance, 1 set of coordinates
      gcdist1(lon = x[,1], lat = x[,2], eps = .Machine$double.eps)
    }
  } else { # if two sets of coordinates
    if (!longlat) {
      # euclidean distance, 2 sets of coordinates
      eucdist2(x1 = x[,1], y1 = x[,2],
               x2 = y[,1], y2 = y[,2],
               eps = .Machine$double.eps)
    } else {
      # great circle distance, 2 sets of coordinates
      gcdist2(lon1 = x[,1], lat1 = x[,2],
              lon2 = y[,1], lat2 = y[,2],
              eps = .Machine$double.eps)
    }
  }
}

#' @rdname gedist
#' @export
eucdist = function(x, y = NULL) {
  arg_check_coords(x, "x")
  if (!is.null(y)) { # if there are 2 sets of coordinates
    arg_check_coords(y, "y")
    eucdist2(x1 = x[,1], y1 = x[,2],
             x2 = y[,1], y2 = y[,2],
             eps = .Machine$double.eps)
  } else { # if there is 1 set of coordinates
    eucdist1(x = x[,1], y = x[,2],
             eps = .Machine$double.eps)
  }
}

#' @rdname gedist
#' @export
gcdist = function(x, y = NULL) {
  arg_check_coords(x, "x")
  if (!is.null(y)) { # if there are 2 sets of coordinates
    arg_check_coords(y, "y")
    gcdist2(lon1 = x[,1], lat1 = x[,2],
            lon2 = y[,1], lat2 = y[,2],
            eps = .Machine$double.eps)

  } else { # if there is 1 set of coordinates
    gcdist1(lon = x[,1], lat = x[,2], eps = .Machine$double.eps)
  }
}

#' Check arguments of a coordinate
#'
#' @param x A 2-d matrix of coordinates
#' @param name The name of the object
#'
#' @noRd
arg_check_coords = function(x, name) {
  if (is.null(dim(x))) {
    stop(paste(name, "must be matrix-like (have non-null dimensions)"))
  }
  dim_x = dim(x)
  if (length(dim_x) != 2) {
    stop(paste(name, "must be matrix-like (have only 2 dimensions)"))
  }
  if (dim_x[1] < 1) {
    stop(paste(name, "must have at least 1 row"))
  }
}
