#' Optimal Population Upper Bound Statistics
#'
#' \code{optimal_ubpop} computes statistics for choosing an optimal population
#' upper bound. \code{ubpop_seq} is a sequence of values to consider as the
#' optimal choice of upper bound. The smallest value must be at least
#' \code{min(pop)/sum(pop)} and should generally be less than or equal to 0.5.
#'
#' @inheritParams scan.test
#' @param ubpop_seq A strictly increasing numeric vector with values between
#'   min(pop)/sum(pop) and 1. The default is \code{seq(0.01, 0.5, len = 50)}.
#'
#' @return Returns a \code{smerc_optimal_ubpop} object. This includes:
#'   \item{ubpop_seq}{The sequence of population bounds considered}
#'   \item{elbow_method}{An object with statistics related to the elbow method}
#'   \item{gini_method}{An object with statistics related to the gini method}
#'   \item{elbow_ubpop}{The population upperbound suggested by the elbow method}
#'   \item{gini_ubpop}{The population upperbound suggested by the Gini method}
#' @seealso \code{\link{scan.test}}
#' @author Joshua French
#' @export
#' @references
#' Meysami, Mohammad, French, Joshua P., and Lipner, Ettie M. The estimation of
#' the optimal cluster upper bound for scan methods in retrospective disease
#' surveillance. Submitted.
#'
#' Han, J., Zhu, L., Kulldorff, M. et al. Using Gini coefficient to determining
#' optimal cluster reporting sizes for spatial scan statistics. Int J Health
#' Geogr 15, 27 (2016). <doi:10.1186/s12942-016-0056-6>
#' @examples
#' data(nydf)
#' coords <- with(nydf, cbind(longitude, latitude))
#' ubpop_stats <- optimal_ubpop(
#'   coords = coords, cases = nydf$cases,
#'   pop = nydf$pop, nsim = 49,
#'   ubpop_seq = seq(0.05, 0.5, by = 0.05)
#' )
#' ubpop_stats
#' \dontrun{
#' plot(ubpop_stats)
#' }
optimal_ubpop <- function(coords, cases, pop,
                          ex = sum(cases) / sum(pop) * pop,
                          nsim = 499, alpha = 0.05,
                          ubpop_seq = seq(0.01, 0.5, len = 50),
                          longlat = FALSE, cl = NULL,
                          type = "poisson",
                          min.cases = 0,
                          simdist = "multinomial") {
  # argument checking
  type <- match.arg(type, c("poisson", "binomial"))
  simdist <- match.arg(simdist, c("multinomial", "poisson", "binomial"))
  arg_check_optimal_ubpop(
    coords = coords, cases = cases,
    pop = pop, ex = ex, nsim = nsim,
    alpha = alpha, ubpop_seq = ubpop_seq,
    longlat = longlat,
    k = 1, w = diag(nrow(coords)),
    type = type, simdist = simdist,
    min.cases = min.cases
  )

  pruned_seq <- seq_scan_test(
    coords = coords,
    cases = cases, pop = pop,
    ex = ex,
    nsim = nsim, alpha = alpha,
    ubpop_seq = ubpop_seq,
    longlat = longlat, cl = cl,
    type = type,
    min.cases = min.cases,
    simdist = simdist
  )
  # compute statistics for elbow method
  estats <- elbow_stats(pruned_seq, ubpop_seq = ubpop_seq)

  # compute statistics for gini method
  gstats <- gini_stats(pruned_seq, ubpop_seq = ubpop_seq, cases = cases, ex = ex, alpha = alpha)

  structure(list(
    ubpop_seq = ubpop_seq,
    elbow_method = estats,
    gini_method = gstats,
    elbow_ubpop = estats$elbow_x,
    gini_ubpop = gstats$gini_x
  ),
  class = "smerc_optimal_ubpop"
  )
}

# compute gini coefficient based on cases in each cluster
# expected cases in each cluster
# total number of cases
gini_coeff <- function(casein, exin, ty) {
  # order by ex
  o <- order(exin, decreasing = FALSE)

  # compute cumulative proportion of cases, ex in sig clusters
  # in order of ex
  cp_casein <- cumsum(casein[o]) / ty
  cp_exin <- cumsum(exin[o]) / ty

  # compute gini coefficient
  # 2 * (.5 - MESS::auc(c(0, cp_casein, 1), c(0, cp_exin, 1)))
  # ensure x and y have appropriate structure
  lin_approx <-
    stats::approx(x = c(0, cp_casein, 1),
                  y = c(0, cp_exin, 1),
                  xout = sort(unique(c(0, cp_casein, 1))))
  # change in x
  dx <- diff(lin_approx$x)
  # sum of y endpoints
  sy <- lin_approx$y[-length(lin_approx$y)] + lin_approx$y[-1]
  # trapezoid rule area
  auc <- sum(dx * sy)/2
  # gini coefficient
  2 * (.5 - auc)
}

#' Compute statistics for elbow method
#'
#' @param x A list returned by the seq_scan_test function
#' @param ubpop_seq A sequence of ubpop values
#'
#' @return An elbow_stats object with the ubpop_seq, stats ( = neg_sum_lrt), elbow_idx, elbow_x, and elbow_y
#' @noRd
elbow_stats <- function(x, ubpop_seq) {
  if (!all.equal(names(x[[1]]), c("tobs", "zones", "pvalue", "idx"))) {
    stop("x must be an object returned by seq_scan_test")
  }
  # extract the significant (or most likely cluster) tobs
  # for each population upper bound
  sig_tobs_seq <- lget(x, name = "tobs")

  # compute the negative sum of the test statistics (or MLC)
  neg_sum_lrt <- -sapply(sig_tobs_seq, sum)

  # compute elbow point
  ep <- elbow_point(ubpop_seq, neg_sum_lrt)

  structure(list(
    ubpop_seq = ubpop_seq,
    stats = neg_sum_lrt,
    elbow_idx = ep$idx,
    elbow_x = ep$x,
    elbow_y = ep$y
  ), class = "elbow_stats")
}

#' Compute gini statistics
#'
#' @param x list returned by seq_scan_test function
#' @param ubpop_seq Sequence of ubpop values
#' @param cases A vector with the expected cases
#' @param ex A vector with the expected counts
#' @param alpha Significance level
#'
#' @return A gini_stats object with the relevant statistics for the gini method
#' @noRd
gini_stats <- function(x, ubpop_seq, cases, ex, alpha) {
  # sanity check x
  if (!all.equal(names(x[[1]]), c("tobs", "zones", "pvalue", "idx"))) {
    stop("x must be an object returned by seq_scan_test")
  }

  # for the significant clusters in x, compute the total cases and total ex
  # get sig pvalues for each population upper bound
  sig_pvalue_seq <- lget(x, name = "pvalue")
  # get sig zones for each population upper bound
  sig_zones_seq <- lget(x, name = "zones")

  # get number of cases in each sig zone for each population upper bound
  sig_yin_seq <- lapply(sig_zones_seq, zones.sum, y = cases)
  # get expected cases in each sig zone for each population upper bound
  sig_ein_seq <- lapply(sig_zones_seq, zones.sum, y = ex)

  # compute gini coefficients for sequence of population upper bounds
  gini_seq <- mapply(gini_coeff,
    casein = sig_yin_seq,
    exin = sig_ein_seq,
    MoreArgs = list(ty = sum(cases))
  )

  # assign any non-significant sums (in case only a MLC returned)
  all_sig <- (sapply(sig_pvalue_seq, min) <= alpha)

  # compute gini coefficient seq corrected
  gini_seq_sig <- gini_seq * all_sig

  # determine index of maximum gini coefficient (first)
  gidx <- which.max(gini_seq)

  structure(list(
    ubpop_seq = ubpop_seq,
    stats = gini_seq_sig,
    gini_idx = gidx,
    gini_x = ubpop_seq[gidx],
    gini_y = gini_seq[gidx]
  ), class = "gini_stats")
}

#' Argument checking for scan tests
#'
#' @param coords A matrix of coordinates
#' @param cases A vector of numeric cases
#' @param pop A vector of population values
#' @param ex A vector of expected counts
#' @param nsim A non-negative integer
#' @param alpha A value greater than 0
#' @param nreport Not used
#' @param ubpop_seq An strictly increasing sequence of values between min(pop)/sum(pop) and 1.
#' @param longlat A logical. TRUE is great circle distance.
#' @param parallel Not used.
#' @param k Number of nearest neighbors. Not always needed.
#' @param w A spatial proximity matrix
#' @param type Statistic type
#' @param simdist Distribution of simulation
#' @param min.cases Minimum number of cases. Only for scan.test.
#' @return NULL
#' @noRd
arg_check_optimal_ubpop <- function(coords, cases, pop, ex, nsim, alpha,
                                    ubpop_seq, longlat, k, w, type, simdist,
                                    min.cases) {
  arg_check_scan_test(
    coords = coords, cases = cases,
    pop = pop, ex = ex, nsim = nsim,
    alpha = alpha, ubpop = 0.1,
    longlat = longlat,
    k = 1, w = diag(nrow(coords)),
    type = type, simdist = simdist,
    min.cases = min.cases
  )
  lb <- min(pop) / sum(pop)
  arg_check_ubpop_seq(ubpop_seq, lb)
}
