#' Check uls.zones arguments
#'
#' @param cases A vector of cases
#' @param pop A vector of populations for each region
#' @param w A spatial adjacency matrix
#' @param ubpop A population upperbound
#' @param check.unique A logical value. TRUE means check for
#' unique rates between regions.
#' @return NULL
#' @noRd
arg_check_uls_zones = function(cases, pop, w, ubpop,
                               check.unique) {
  N = length(cases)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_w(w, N)
  arg_check_ubpop(ubpop)
  arg_check_check_unique(check.unique)
}
