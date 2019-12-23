#' Check arguments for \code{smerc_cluster}
#'
#' \code{smerc_cluster} prepares and returns a
#' \code{smerc_cluster}.
#'
#' @param tobs The vector of observed test statistics for each zone
#' @param zones A list of zones
#' @param pvalue The p-value associated with each test statistic
#' @inheritParams flex.test
#' @param d A precomputed distance matrix based on \code{coords}
#' @param method A character string indicating the method
#' used to construct the \code{smerc_cluster}.
#' @param rel_param A names list with the relevant parameters
#' associated with \code{method}.
#' @param a A single value >= 0 indicating the penalty to use
#' for \code{\link{elliptic.test}}.
#' @param shape_all A vector of shape parameters associated
#' with \code{zones}. Relevant for \code{\link{elliptic.test}}.
#' @param angle_all A vector of angle parameter associated with
#' \code{zones}. Relevant for \code{\link{elliptic.test}}.
#' @param alpha The significance level.
#' @return A \code{smerc_cluster} object
#' @noRd
arg_check_smerc_cluster = function(tobs, zones, pvalue,
                                   coords, cases, pop, ex,
                                   longlat, method,
                                   rel_param, w, d, a,
                                   shape_all, angle_all,
                                   alpha) {
  arg_check_tobs(tobs)
  Ntobs = length(tobs)
  arg_check_zones(zones, Ntobs)
  arg_check_pvalue(pvalue, Ntobs)
  arg_check_coords(coords)
  N = nrow(coords)
  arg_check_cases(cases, N)
  arg_check_pop(pop, N)
  arg_check_ex(ex, N)
  arg_check_longlat(longlat)
  arg_check_method(method)
  arg_check_rel_param(rel_param)
  if (!is.null(w)) {
    arg_check_w(w, N)
  }
  if (!is.null(d)) {
    arg_check_d(d, N)
  }
  if (!is.null(a)) {
    arg_check_a(a)
  }
  if (!is.null(shape_all)) {
    arg_check_shape_all(shape_all, Ntobs)
  }
  if (!is.null(angle_all)) {
    arg_check_angle_all(angle_all, Ntobs)
  }
  arg_check_alpha(alpha)
}
