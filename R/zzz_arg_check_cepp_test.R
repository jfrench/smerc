#' Argument checking for cepp.test
#'
#' Check the arguments of the cepp.test function
#' @return NULL
#' @noRD
arg_check_cepp_test = function(coords, cases, pop, nstar,
                               ex, nsim, longlat, alpha, noc,
                               simdist) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  if (length(nstar) != 1 || !is.numeric(nstar)) {
    stop("nstar should be a numeric vector of length 1")
  }
  if (nstar < 1 || nstar > sum(pop)) {
    stop("nstar should be at least 1 and less than or equal to the sum(pop)")
  }
  arg_check_ex(ex, N)
  arg_check_alpha(alpha)
  arg_check_nsim(nsim)
  arg_check_longlat(longlat)
  arg_check_noc(noc)
  arg_check_simdist(simdist)
}
