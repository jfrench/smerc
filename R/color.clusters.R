#' Color clusters
#'
#' \code{color.clusters} is a helper function to color
#' clusters of regions produced by an appropriate method,
#' e.g., \code{scan.test} or \code{uls.test}.  Regions that
#' are not part of any cluster have no color.
#'
#' @param x An object of class scan produced by a function
#'   such as \code{scan.test}.
#' @param col A vector of colors to color the clusters in
#'   \code{x}.  Should have same length as the number of
#'   clusters in \code{x}.
#' @inheritParams summary.smerc_cluster
#' @return Returns a vector with colors for each
#'   region/centroid for the data set used to construct
#'   \code{x}.
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' coords <- with(nydf, cbind(longitude, latitude))
#' out <- scan.test(
#'   coords = coords, cases = floor(nydf$cases),
#'   pop = nydf$pop, alpha = 0.2, longlat = TRUE,
#'   nsim = 9
#' )
#' data(nypoly)
#' library(sp)
#' # plot all clusters
#' plot(nypoly, col = color.clusters(out), axes = TRUE)
#' # zoom in on small cluster
#' plot(nypoly,
#'   col = color.clusters(out),
#'   xlim = c(400000, 450000),
#'   ylim = c(4750000, 4800000)
#' )
#' # plot only clusters 1 and 3
#' plot(nypoly, col = color.clusters(out, idx = c(1, 3)))
color.clusters <- function(x, idx = seq_along(x$clusters), col = grDevices::hcl.colors(length(idx))) {
  if (!any(is.element(class(x), c("scan", "smerc_cluster")))) {
    stop("x should be an object of class scan or smerc_cluster.")
  }
  if (min(idx) < 1 | max(idx) > length(x$clusters)) {
    stop("invalid idx values")
  }
  if (length(idx) != length(col)) {
    stop("The number of colors must match the length of idx.")
  }

  mycol <- numeric(nrow(x$coords))
  for (i in seq_along(x$clusters)) {
    mycol[x$clusters[[i]]$loc] <- col[i]
  }
  return(mycol)
}
