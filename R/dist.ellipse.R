#' Compute minor axis distance of ellipse
#'
#' \code{dist.ellipse} computes the length of the minor axis needed for an
#' ellipse of a certain \code{shape} and \code{angle} to intersect
#' each of the other coordinates from a starting coordinate.
#'
#' @param coords An \eqn{N \times 2} matrix of coordinates
#' @param shape The ratio of the major axis to the minor axis of the ellipse
#' @param angle The angle of the ellipse in the range [0, 180).
#'
#' @return A matrix of distances between each coordinate and all other coordinates (and itself).  Each row contains the distances for a coordinate.
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' d = dist.ellipse(coords, 4, 15)
dist.ellipse = function(coords, shape, angle) {
  arg_check_dist_ellipse(coords, shape, angle)
  rad = angle * pi/180
  N = nrow(coords)
  m3 = matrix(c(sin(rad), cos(rad)), ncol = 2, 
              nrow = N, byrow = TRUE)
  t(sapply(seq_len(N), function(i) {
    coordsi = matrix(coords[i,], byrow = TRUE, 
                     ncol = 2, nrow = N)
    cdiff = (coords - coordsi)
    sqrt(((cdiff[,1] * m3[,2] + cdiff[,2] * m3[,1])^2 + 
           (cdiff[,1] * m3[,1] - cdiff[,2] * m3[,2])^2 * shape^2)/shape^2)
    
  }))
}

arg_check_dist_ellipse = function(coords, shape, angle) {
  if (!is.matrix(coords)) {
    stop("coords must be a matrix")
  }
  if (length(shape) != 1) {
    stop("shape should be a single value")
  }
  if (shape < 1) {
    stop("shape should be at least 1")
  }
  if (length(angle) != 1) {
    stop("angle should be a single value")
  }
  if (angle < 00 | angle >= 360) {
    stop("angle should be in [0, 360) degrees")
  }
}

all.shape.dists = function(s, na, coords) {
  angle = seq(90, 270, len = na + 1)[-(na + 1)]
  d = lapply(seq_along(angle), function(j) {
    dist.ellipse(coords, shape = s, angle = angle[j])
  })
  do.call(rbind, d)
}

