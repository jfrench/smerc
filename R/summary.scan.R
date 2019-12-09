#' Summary of \code{scan} object
#'
#' Summary of object of class \code{scan}.
#'
#' @param x An object of class \code{scan}.
#' @inheritDotParams base::summary
#' @param idx A index vector indicating the elements of
#' \code{x$clusters} to print information for. The default
#' is all clusters.
#' @param digits Integer indicating the number of decimal places.
#' @method summary scan
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' summary(out)
summary.scan = function(x, ...,
                        idx = seq_along(x$clusters),
                        digits = 1) {
  if (class(x) != "scan") {
    stop("x must be a scan object.")
  }
  if (min(idx) < 1 | max(idx) > length(x$clusters)) {
    stop("invalid idx values")
  }
  regions = sapply(x$clusters[idx], function(i) length(i$locids))
  max_dist = base::round(sget(x$clusters, "max_dist"),
                         digits = digits)
  cases = sget(x$clusters[idx], "cases")
  ex = base::round(sget(x$clusters[idx], "expected"),
                   digits = digits)
  rr = base::round(sget(x$clusters[idx], "rr"),
                   digits = digits)
  stat = base::round(sget(x$clusters[idx], "test_statistic"),
                     digits = 1)
  p = sget(x$clusters[idx], "pvalue")
  data.frame(regions,
             max_dist,
             cases,
             ex,
             rr,
             stat,
             p
  )
}

