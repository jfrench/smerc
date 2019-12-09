#' Print object of class \code{scan}.
#'
#' Print \code{scan} object
#'
#' @param x An object of class \code{scan}.
#' @inheritDotParams base::print
#' @method print scan
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' out
print.scan = function(x, ...) {
  message(paste("Total regions:", x$number_of_regions))
  message(paste("Total cases:", x$total_cases))
  message(paste("Total population:", x$total_population))
  message(paste("Cases per 100,000 persons:", x$cases_per_100k))
}

