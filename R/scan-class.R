#' Class\code{scan}.
#'
#' An object of class \code{scan}.  A \code{scan} object
#' generally has the following components:
#' \item{clusters}{A list containing information about the
#' significant clusters.}
#' \item{coords}{A matrix with the centroid coordinates of
#' the observed regions.}
#' \item{number_of_regions}{The number of observed regions.}
#' \item{total_population}{The total population across all regions.}
#' \item{total_cases}{The total number of cases across all regions.}
#' \item{cases_per_100k}{The estimated number of cases for a population of size 100,000.}
#' The \code{clusters} component generally has the following elements, though some
#' may be \code{NA} or missing if not appropriate:
#' \item{locids}{The indices of the regions in the significant cluster.}
#' \item{coords}{The originating centroid of the significant cluster.}
#' \item{r}{The radius of the window of the clusters.}
#' \item{max_dist}{The maximum distance between the centroids of the significant cluster.}
#' \item{pop}{The total population in the cluster.}
#' \item{cases}{The observed number of cases in the cluster.}
#' \item{expected}{The expected number of cases in the cluster.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{loglikrat}{The log of likelihood ratio test statistic for the cluster.}
#' \item{test_statistic}{The test statistic for the cluster.}
#' \item{pvalue}{The p-value of the test statistic associated with the cluster window.}
#' \item{w}{The connectivity information for the cluster}.
#' For \code{\link{elliptic.test}}, \code{cluster} additionally has:
#' \item{semiminor_axis}{The semi-minor axis length for the ellipse.}
#' \item{semimajor_axis}{The semi-major axis length for the ellipse.}
#' \item{angle}{The rotation angle of the ellipse.}
#' \item{shape}{The shape of the ellipse.}
#' @name scan_class
#' @rdname scan_class
#' @exportClass scan
#' @export
#'
