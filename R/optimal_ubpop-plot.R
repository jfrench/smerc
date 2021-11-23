#' Plot object of class \code{smerc_optimal_ubpop}.
#'
#' Plot results of \code{\link{optimal_ubpop}}. This is only meant for a visual
#' summary of the results. Users will need to access the elements of the \code{smerc_optimal_ubpop} object \code{x}
#' if they want to create a custom plot.
#'
#' @param x An object of class \code{smerc_optimal_ubpop}.
#' @param ... Not used
#' @param method The method to plot. The default is \code{"all"}.
#' The other valid options are \code{"elbow"} and \code{"gini"}.
#' @export
#' @seealso \code{\link{optimal_ubpop}}
#' @method plot smerc_optimal_ubpop
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' ubpop_stats = optimal_ubpop(coords = coords, cases  = nydf$cases,
#'                             pop = nydf$pop, nsim = 49,
#'                             ubpop = seq(0.05, 0.5, by = 0.05))
#' \dontrun{plot(ubpop_stats)}
#' plot(ubpop_stats, method = "elbow")
#' plot(ubpop_stats$ubpop_seq, ubpop_stats$elbow_method$stats)
#' plot(ubpop_stats, method = "gini")
#' plot(ubpop_stats$ubpop_seq, ubpop_stats$gini_method$stats)
plot.smerc_optimal_ubpop = function(x, ..., method = "all") {
  if (method == "all") {
    # get current value of mfrow
    original_mfrow = graphics::par()$mfrow
    graphics::par(mfrow = c(1, 2))
    plot_elbow(x)
    plot_gini(x)
    # return mfrow to original mfrow
    graphics::par(mfrow = original_mfrow)
  } else if (method == "elbow") {
    plot_elbow(x)
  } else if (method == "gini") {
    plot_gini(x)
  }
}

plot_elbow = function(x) {
  graphics::plot(x$ubpop_seq, x$elbow_method$stats, xlab = "% population", ylab = "-LRT")
  graphics::points(x$elbow_method$elbow_x, x$elbow_method$elbow_y, pch = 19)
  graphics::abline(v = x$elbow_ubpop)
  graphics::title(paste("elbow:", round(x$elbow_ubpop, 2) * 100, "%"))
}

plot_gini = function(x) {
  graphics::plot(x$ubpop_seq, x$gini_method$stats, xlab = "% population", ylab = "gini coefficient")
  graphics::abline(v = x$gini_ubpop)
  graphics::title(paste("gini:", round(x$gini_ubpop, 2) * 100, "%"))
}
