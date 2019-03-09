#' Argument checking for many scan-like functions
#'
#' Check the arguments of many \code{*.test} functions.
#' @return NULL
#' @export
#' @keywords internal
# argument checking for most scan tests
arg_check_scan_test =
  function(coords, cases, pop, ex, nsim, alpha, nreport,
           ubpop, longlat, parallel, k, w, type = NULL) {
    if (!(is.matrix(coords) | is.data.frame(coords))) {
      stop("coords should be a matrix or a data frame")
    }
    if (ncol(coords) != 2) {
      stop("coords must have two columns")
    }
    N = nrow(coords)
    if (length(cases) != N) {
      stop("length(cases) != nrow(coords)")
    }
    if (!is.numeric(cases)) {
      stop("cases should be a numeric vector")
    }
    if (length(pop) != N) {
      stop("length(pop) != nrow(coords)")
    }
    if (!is.numeric(pop)) {
      stop("pop should be a numeric vector")
    }
    if (length(ex) != N) {
      stop("length(ex) != nrow(coords)")
    }
    if (!is.numeric(ex)) {
      stop("ex should be a numeric vector")
    }
    if (length(alpha) != 1 || !is.numeric(alpha)) {
      stop("alpha should be a numeric vector of length 1")
    }
    if (alpha < 0 || alpha > 1) {
      stop("alpha should be a value between 0 and 1")
    }
    if (length(nsim) != 1 || !is.numeric(nsim)) {
      stop("nsim should be a vector of length 1")
    }
    if (nsim < 0) {
      stop("nsim should be an non-negative integer")
    }
    if (length(ubpop) != 1 || !is.numeric(ubpop)) {
      stop("ubpop should be a numeric vector of length 1")
    }
    if (ubpop <= 0 || ubpop > 1) {
      stop("ubpop should be a value between 0 and 1")
    }
    if (length(longlat) != 1) {
      stop("length(longlat) != 1")
    }
    if (!is.logical(longlat)) {
      stop("longlat should be a logical value")
    }
    if (length(parallel) != 1) {
      stop("length(parallel) != 1")
    }
    if (!is.logical(parallel)) {
      stop("parallel should be a logical value")
    }
    if (length(k) != 1) {
      stop("k must have length 1")
    }
    if (k < 1) stop("k must be an integer >= 1")
    if (!is.matrix(w)) {
      stop("w must be a matrix")
    }
    if (nrow(w) != ncol(w)) {
      stop("w much be a square matrix")
    }
    if (!is.numeric(w)) {
      stop("w must be a numeric matrix")
    }
    if (nrow(w) != nrow(coords)) {
      stop("nrow(w) != nrow(coords)")
    }
    if (floor(k) > nrow(coords)) {
      stop("k cannot be more than the number of regions")
    }
    if (!is.null(type)) {
      if (!is.element(type, c("poisson", "binomial"))) {
        stop("type must be 'poisson' or 'binomial'")
      }
    }
  }
