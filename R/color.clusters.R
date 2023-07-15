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
#' set.seed(1)
#' data(nydf)
#' coords <- with(nydf, cbind(longitude, latitude))
#' out <- scan.test(
#'   coords = coords, cases = floor(nydf$cases),
#'   pop = nydf$pop, alpha = 0.2, longlat = TRUE,
#'   nsim = 9
#' )
#' #' # better plotting
#' if (require("sf", quietly = TRUE)) {
#'    data(nysf)
#'    plot(st_geometry(nysf), col = color.clusters(out))
#'    # plot only clusters 2 and 3
#'    plot(st_geometry(nysf),
#'         col = color.clusters(out, idx = c(2, 3)),
#'         border = "white")
#' }
color.clusters <- function(x, idx = seq_along(x$clusters),
                           col = grDevices::hcl.colors(length(idx))) {
  if (!any(is.element(class(x), c("scan", "smerc_cluster")))) {
    stop("x should be an object of class scan or smerc_cluster.")
  }
  if (min(idx) < 1 | max(idx) > length(x$clusters)) {
    stop("invalid idx values")
  }
  if (length(idx) != length(col)) {
    stop("The number of colors must match the length of idx.")
  }

  mycol <- rep(NA, nrow(x$coords))
  for (i in seq_along(idx)) {
    mycol[x$clusters[[idx[i]]]$loc] <- col[i]
  }
  return(mycol)
}
