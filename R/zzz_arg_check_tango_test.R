#' Argument checking for tango.test
#'
#' Check the arguments of the tango.test function
#' @return NULL
#' @noRd
arg_check_tango_test = function(cases, pop, w, nsim) {
  N = length(cases)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_tango_w(w, N)
  arg_check_nsim(nsim)
}
