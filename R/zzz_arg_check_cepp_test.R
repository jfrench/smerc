#' Argument checking for cepp.test
#'
#' Check the arguments of the cepp.test function
#'
#' @param coords A matrix of coordinates
#' @param cases A vector of case counts
#' @param pop A vector of population values
#' @param nstar Window radius (in terms of population)
#' @param ex A vector of expected counts
#' @param nsim Number of simulations to perform
#' @param longlat Logical. TRUE = great circle distance
#' @param alpha Significance level
#' @param simdist Simulation distribution
#'
#' @return NULL
#' @noRd
arg_check_cepp_test = function(coords, cases, pop, nstar,
                               ex, nsim, longlat, alpha,
                               simdist) {
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_nstar(nstar, pop)
  arg_check_ex(ex, N)
  arg_check_alpha(alpha)
  arg_check_nsim(nsim)
  arg_check_longlat(longlat)
  arg_check_simdist(simdist)
}
