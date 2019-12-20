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
  arg_check_nstar(nstar)
  arg_check_ex(ex, N)
  arg_check_alpha(alpha)
  arg_check_nsim(nsim)
  arg_check_longlat(longlat)
  arg_check_noc(noc)
  arg_check_simdist(simdist)
}
