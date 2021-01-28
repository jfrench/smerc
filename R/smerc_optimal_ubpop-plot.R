#' Plot object of class \code{smerc_optimal_ubpop}.
#'
#' Plot results of \link{\code{optimal_ubpop}}. This is only meant for a visual
#' summary of the results. Users will need to access the elements of the \code{smerc_optimal_ubpop} object \code{x}
#' if they want to create a custom plot.
#'
#' @param x An object of class \code{smerc_optimal_ubpop}.
#' @param ... Not used
#' @param method The method to plot. The default is \code{"all"}.
#' The other valid options are \code{"elbow"} and \code{"gini"}.
#' @export
#' @method plot smerc_optimal_ubpop
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' ubpop_stats = optimal_ubpop(coords = coords, cases  = nydf$cases,
#'                             pop = nydf$pop, nsim = 49,
#'                             ubpop = seq(0.05, 0.5, by = 0.05))
#' plot(ubpop_stats)
#' plot(ubpop_stats, method = "elbow")
#' plot(ubpop_stats$ubpop_seq, ubpop_stats$neg_lrt)
#' plot(ubpop_stats, method = "gini")
#' plot(ubpop_stats$ubpop_seq, ubpop_stats$gini_coef)
plot.smerc_optimal_ubpop = function(x, ..., method = "all") {
  if (method == "all") {
    graphics::par(mfrow = c(1, 2))
    plot_elbow(x)
    plot_gini(x)
  } else if (method == "elbow") {
    plot_elbow(x)
  } else if (method == "gini") {
    plot_gini(x)
  }
}

plot_elbow = function(x) {
  graphics::plot(x$ubpop_seq, x$neg_lrt, xlab = "% population", ylab = "-LRT")
  graphics::points(x$eb_point$x, x$eb_point$y, pch = 19)
  graphics::abline(v = x$elbow_ubpop)
  graphics::title(paste("elbow:", round(x$elbow_ubpop, 2) * 100, "%"))
}

plot_gini = function(x) {
  graphics::plot(x$ubpop_seq, x$gini_coef, xlab = "% population", ylab = "gini coefficient")
  graphics::abline(v = x$gini_ubpop)
  graphics::title(paste("gini:", round(x$gini_ubpop, 2) * 100, "%"))
}
