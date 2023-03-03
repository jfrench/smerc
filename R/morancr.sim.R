#' Constant-risk Moran's I statistic
#'
#' \code{morancr.stat} computes the constant-risk version of the Moran's I
#' statistic proposed by Walter (1992).
#'
#' @inheritParams flex.test
#'
#' @return Returns a numeric value.
#' @author Joshua French
#' @seealso \code{\link{morancr.test}}
#' @export
#' @references Walter, S. D. (1992). The analysis of regional patterns in health
#'   data: I. Distributional considerations. American Journal of Epidemiology,
#'   136(6), 730-741.
#' @examples
#' data(nydf)
#' data(nyw)
#' ex <- sum(nydf$cases) / sum(nydf$pop) * nydf$pop
#' morancr.sim(nsim = 10, cases = nydf$cases, w = nyw, ex = ex)
morancr.sim <- function(nsim = 1, cases, w, ex) {
  simdist <- "multinomial"
  ty <- sum(cases)
  # compute max test stat for nsim simulated data sets
  tsim <- pbapply::pblapply(seq_len(nsim), function(i) {
    # simulate new data
    if (simdist == "multinomial") {
      ysim <- stats::rmultinom(1, size = ty, prob = ex)
    }
    y_std <- matrix((ysim - ex) / sqrt(ex))
    return(sum(w * y_std %*% t(y_std)) / sum(w))
  })
  unlist(tsim, use.names = FALSE)
}
