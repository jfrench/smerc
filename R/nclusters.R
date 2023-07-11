#' Number of clusters
#'
#' \code{nclusters} returns the number of clusters
#' identified in a \code{smerc_cluster} object.
#' @param x A \code{smerc_cluster object}
#'
#' @return A non-negative integer.
#' @export
#' @examples
#' data(nydf)
#' coords <- with(nydf, cbind(longitude, latitude))
#' out <- scan.test(
#'   coords = coords, cases = floor(nydf$cases),
#'   pop = nydf$pop, nsim = 19,
#'   alpha = 0.3, longlat = TRUE
#' )
#' nclusters(out)
nclusters <- function(x) {
  if (!is.element("smerc_cluster", class(x))) {
    stop("x must be a smerc_cluster")
  }
  length(x$clusters)
}
