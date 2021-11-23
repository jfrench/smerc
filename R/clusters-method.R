#' Extract clusters
#'
#' \code{clusters} extracts the clusters contained in \code{x}.
#'
#' @param x An object with clusters.
#' @param idx An index vector indicating the elements of
#' \code{x$clusters} to print information for. The default
#' is all clusters.
#' @param ... Currently unimplemented
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 19,
#'                 alpha = 0.2, longlat = TRUE)
#' clusters(out)
#' clusters(out, idx = 1:2)
clusters = function(x, idx = seq_along(x$clusters), ...) {
  UseMethod("clusters")
}
