#' Extract clusters
#'
#' \code{clusters} extracts the clusters contained in a
#' \code{\link{smerc_cluster}} object.
#'
#' @param x An object of class \code{smerc_cluster}.
#'
#' @return A list. Each element of the list is a vector with the indices of event
#' locations in the associated cluster.
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 19,
#'                 alpha = 0.2, longlat = TRUE)
#' clusters(out)
clusters = function(x) {
  if (!is.element("smerc_cluster", class(x))) {
    stop("x must have smerc_cluster class")
  }
  lget(x$clusters, "locids")
}
