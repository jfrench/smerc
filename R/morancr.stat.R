#' Constant-risk Moran's I statistic
#'
#' \code{morancr.stat} computes the constant-risk version of the Moran's I
#' statistic proposed by Walter (1992).
#'
#' @param ex The expected number of cases for each region.
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
#' morancr.stat(cases = nydf$cases, w = nyw, ex = ex)
morancr.stat <- function(cases, w, ex) {
  arg_check_morancr_stat(cases = cases, pop = cases + 1, w = w, ex = ex)
  y_std <- matrix((cases - ex) / sqrt(ex))
  return(sum(w * y_std %*% t(y_std)) / sum(w))
}

#' Argument checking for moran_cr.test
#'
#' @noRd
arg_check_morancr_stat <- function(cases, pop, w, ex) {
  if (!is.numeric(cases)) {
    stop("cases should be a numeric values")
  }
  if (!is.vector(cases)) {
    stop("cases should be a vector")
  }
  if (min(cases) < 0) {
    stop("cases must have non-negative values")
  }
  N <- length(cases)
  if (length(pop) != N) {
    stop("length(pop) != length(cases)")
  }
  if (!is.numeric(pop)) {
    stop("pop should be numeric values")
  }
  if (!is.vector(pop)) {
    stop("pop should be a vector")
  }
  if (min(pop) < 0) {
    stop("pop values must be >= 0")
  }
  if (length(ex) != N) {
    stop("length(ex) != length(cases)")
  }
  if (!is.numeric(ex)) {
    stop("ex should be numeric values")
  }
  if (!is.vector(ex)) {
    stop("ex should be a vector")
  }
  if (!is.matrix(w) & !is.data.frame(w)) {
    stop("w must be a matrix or data.frame")
  }
  if (nrow(w) != N | ncol(w) != N) {
    stop("w must be a square matrix with nrow(w) = length(coords)")
  }
  if (any(w != 0 & w != 1)) {
    stop("w must be 0s and 1s")
  }
}
