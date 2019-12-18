#' Argument checking for scan tests
#'
#' @param coords A matrix of coordinates
#' @param cases A vector of numeric cases
#' @param pop A vector of population values
#' @param ex A vector of expected counts
#' @param nsim A non-negative integer
#' @param alpha A value greater than 0
#' @param nreport Not used
#' @param ubpop A value between 0 and 1
#' @param longlat A logical. TRUE is great circle distance.
#' @param parallel Not used.
#' @param k Number of nearest neighbors. Not always needed.
#' @param w A spatial proximity matrix
#' @param type Statistic type
#' @param simdist Distribution of simulation
#' @param min.cases Minimum number of cases. Only for scan.test.
#' @return NULL
#' @noRd
arg_check_scan_test =
  function(coords, cases, pop, ex, nsim, alpha,
           nreport = NULL,
           ubpop, longlat, parallel = NULL, k, w, type = NULL,
           simdist = NULL, min.cases = NULL) {
    arg_check_coords(coords)
    N = nrow(coords)
    arg_check_cases(cases, N)
    arg_check_pop(pop, N)
    arg_check_ex(ex, N)
    arg_check_nsim(nsim)
    arg_check_alpha(alpha)
    # nreport no check, deprecated
    arg_check_ubpop(ubpop)
    arg_check_longlat(longlat)
    # parallel no check, deprecated
    arg_check_k(k, N)
    arg_check_w(w, N)
    if (!is.null(type)) {
      arg_check_type(type)
    }
    if (!is.null(simdist)) {
      arg_check_simdist(simdist)
    }
    arg_check_simdist(simdist)
    if (!is.null(min.cases)) {
      arg_check_min_cases(min.cases)
    }
  }
