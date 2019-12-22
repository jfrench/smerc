#' Argument checking for *.sim functions
#'
#' Check the arguments of the \code{*.sim} functions.
#'
#' @param nsim Number of simulations
#' @param ty Total number of cases
#' @param ex Expected counts
#' @param type Type of statistic
#' @param nn List of nn (e.g., nnpop function)
#' @param zones List of zones (e.g., scan.zones)
#' @param ein List of expected in each zone
#' @param eout List of expected out of each zone
#' @param tpop Total population
#' @param popin Population in each zone
#' @param popout Population outside of each zone
#' @param w Spatial adjacency matrix
#' @param pop Vector of populations
#' @param ubpop Population upperbound
#' @param static Static zones. Logical. TRUE for scan.test.
#' FALSE for uls.test.
#' @param simdist Simulation distribution.
#' @return NULL
#' @noRd
arg_check_sim = function(nsim, ty, ex, type,
                         nn = NULL, zones = NULL,
                         ein = NULL, eout = NULL,
                         tpop = NULL, popin = NULL,
                         popout = NULL, w = NULL,
                         pop = NULL, ubpop = NULL,
                         static = FALSE,
                         simdist = "multinomial") {
  arg_check_nsim(nsim)
  arg_check_ty(ty)
  N = length(ex)
  arg_check_ex(ex, N)
  arg_check_type(type)
  if (!is.null(nn)) {
    if (!is.list(nn)) stop("nn must be a list")
    nz = sum(sapply(nn, length))
  }
  if (!is.null(zones)) {
    if (!is.list(zones)) stop("zones must be a list")
    nz = length(zones)
  }
  # if (type == "poisson") {
  #   arg_check_sim_poisson_type(ein = ein, eout = eout, nz = nz)
  # } else if (type == "binomial") {
  #   arg_check_sim_binomial_type(popin = popin, popout = popout,
  #                               tpop = tpop, nz = nz)
  # }
  arg_check_w(w, N)
  if (!is.null(ubpop)) {
    arg_check_ubpop
  }
  if (!is.null(pop)) {
    arg_check_pop(pop, N)
  }
  arg_check_simdist(simdist)
  if (simdist == "binomial" & is.null(pop)) {
    stop("pop must be specified when simdist == 'binomial'")
  }
}
