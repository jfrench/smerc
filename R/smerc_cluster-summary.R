#' Summary of \code{smerc_cluster} object
#'
#' Summary of object of class \code{smerc_cluster}.
#'
#' @param object An object of class \code{smerc_cluster}.
#' @inheritDotParams base::summary
#' @param idx An index vector indicating the elements of
#' \code{object$clusters} to print information for. The default
#' is all clusters.
#' @param digits Integer indicating the number of decimal places.
#' @method summary smerc_cluster
#' @return A \code{data.frame} with columns:
#' \item{nregions}{The number of regions in the cluster.}
#' \item{max_dist}{The maximum intercentroid distance between all the regions
#' in the cluster.}
#' \item{cases}{The number of cases in the cluster.}
#' \item{ex}{The expected number of cases in the cluster.}
#' \item{rr}{Relative risk in the cluster window. This is
#' \code{(cases/pop)/((total_cases - cases)/
#' (total_population - population))}.}
#' \item{stat}{The test statistic for the cluster.}
#' \item{p}{The p-value of the test statistic associated with the cluster.}
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' summary(out)
summary.smerc_cluster = function(object, ...,
                                 idx = seq_along(object$clusters),
                                 digits = 1) {
  if (min(idx) < 1 | max(idx) > length(object$clusters)) {
    stop("invalid idx values")
  }
  nregions = sapply(object$clusters[idx], function(i) length(i$locids))
  max_dist = base::round(sget(object$clusters[idx], "max_dist"),
                         digits = digits)
  cases = sget(object$clusters[idx], "cases")
  ex = base::round(sget(object$clusters[idx], "expected"),
                   digits = digits)
  rr = base::round(sget(object$clusters[idx], "rr"),
                   digits = digits)
  stat = base::round(sget(object$clusters[idx], "test_statistic"),
                     digits = 1)
  p = sget(object$clusters[idx], "pvalue")
  data.frame(nregions,
             max_dist,
             cases,
             ex,
             rr,
             stat,
             p
  )
}
