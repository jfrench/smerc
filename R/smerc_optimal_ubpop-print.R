#' Print object of class \code{smerc_optimal_ubpop}.
#'
#' Print \code{smerc_optimal_ubpop} object
#'
#' @param x An object of class \code{smerc_optimal_ubpop}.
#' @inheritDotParams base::print
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' ubpop_stats = optimal_ubpop(coords = coords, cases  = nydf$cases,
#'                             pop = nydf$pop, nsim = 49,
#'                             ubpop = seq(0.05, 0.5, by = 0.05))
#' ubpop_stats
print.smerc_optimal_ubpop = function(x, ...) {
  cat("ubpop recommendations:\n")
  print(data.frame(method = c("elbow", "gini"),
             ubpop = c(x$elbow_ubpop, x$gini_ubpop)))
  # cat("elbow method recommends ubpop = ", x$elbow_ubpop, "\n")
  # cat("gini method recommends ubpop = ", x$gini_ubpop, "\n")
}
