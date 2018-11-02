#' Spatial scan statistic
#'
#' \code{scan.stat2} calculates the spatial scan statistic
#' for a zone (a set of spatial regions).  The statistic is
#' the log of the likelihood ratio test statistic of the
#' chosen distribution.  If \code{type = "poisson"} and
#' \code{a} is more than zero, this statistic is penalized.
#'
#' @param yin The total number of cases in the zone.
#' @param ty The total number of cases in the study area.
#' @param ein The expected number of cases in the zone.
#'   Conventionally, this is the estimated overall disease
#'   risk across the study area, multiplied by the total
#'   population size of the zone.
#' @param tpop The total population in the study area.
#' @param popin The total population in the zone.
#' @param type The type of scan statistic to implement. The
#'   default choice are \code{"poisson"}.  The other choice
#'   is \code{"binomial"}.
#' @param shape The shape of the ellipse, which is the ratio
#'   of the length of the longest and shortest axes of the
#'   ellipse.  The default is 1, meaning it is a circle.
#' @param a A tuning parameter for the adjusted
#'   log-likelihood ratio.  See details.
#' @param yout The observed number of cases outside the
#'   zone.  This should be \code{ty - yin} and is computed
#'   automatically if not provided.
#' @param eout The expected number of cases outside the
#'   zone.  This should be \code{ty - ein} and is computed
#'   automatically if not provided.
#' @param popout The population outside the zone.  This
#'   should be \code{tpop - popin} and is computed
#'   automatically if not provided.
#' @return A vector of scan statistics.
#' @author Joshua French
#' @export
#' @references
#'
#' Poisson scan statistic:  Kulldorff, M. (1997) A spatial
#' scan statistic. Communications in Statistics - Theory and
#' Methods, 26(6): 1481-1496,
#' <doi:10.1080/03610929708831995>
#'
#' Penalized Poisson scan statistic: Kulldorff, M., Huang,
#' L., Pickle, L. and Duczmal, L. (2006) An elliptic spatial
#' scan statistic. Statistics in Medicine, 25:3929-3943.
#' <doi:10.1002/sim.2490>
#'
#' Binomial scan statistic: Duczmal, L. and Assuncao, R.
#' (2004) A simulated annealing strategy for the detection
#' of arbitrarily shaped spatial clusters. Computational
#' Statistics & Data Analysis, 45(2):269-286.
#' <doi:10.1016/S0167-9473(02)00302-X>
#' @examples
#' # New York leukemia data
#' # total cases
#' ty = 552
#' # total population
#' tpop = 1057673
#'
#' # poisson example with yin = 106 and ein = 62.13
#' scan_stat(yin = 106, ty = ty, ein = 62.13)
#' stat.poisson(yin = 106, yout = 552 - 106,
#'              ein = 62.13, eout = 552 - 62.13)
#'
#' # binomial example with yin = 41 and popin = 38999
#' scan_stat(yin = 41, ty = ty,
#'           popin = 38999, tpop = tpop, type = "binomial")
#' stat.binom(41, ty - 41, ty, 38999, tpop - 38999, tpop)
scan_stat = function(yin, ty, ein = NULL, eout = NULL,
                     popin = NULL, tpop = NULL,
                     type = "poisson", a = 0, shape = 1,
                     yout = NULL,
                     popout = NULL) {
  arg_check_scan_stat(yin = yin, ty = ty, ein = ein, 
                      popin = popin, tpop = tpop, 
                      type = type, a = a, shape = shape)
  B = length(yin)
  # double check yout argument
  if (is.null(yout)) yout = ty - yin
  if (B != length(yout)) stop("length(yin) != length(yout)")
  
  if (type == "poisson") {
    # check eout argument
    if (is.null(eout)) eout = ty - ein 
    if (B != length(eout)) stop("length(yin) != length(eout)")

    tall = stat.poisson(yin, yout, ein, eout, a, shape)
    # tall = yin * (log(yin) - log(ein)) + yout * (log(yout) - log(eout))
    # # correct test statistics for NaNs
    # tall[yin/ein <= yout/eout | is.nan(tall)] = 0
    # if (a > 0) {
    #   wp = which(shape > 1)
    #   tall[wp] = tall[wp] * ((4 * shape[wp])/(shape[wp] + 1)^2)^a
    # }
  } else if (type == "binomial") {
    # check eout argument
    if (is.null(popout)) popout = tpop - popin 
    if (B != length(popout)) stop("length(yin) != length(popout)")
    
    tall = stat.binom(yin, yout, ty, popin, popout, tpop)
  }
  return(tall)
}

#' @export
#' @rdname scan_stat
stat.poisson = function(yin, yout, ein, eout, a = 0, shape = 1) {
  tall = yin * (log(yin) - log(ein)) + yout * (log(yout) - log(eout))
  # correct test statistics for NaNs
  tall[yin/ein <= yout/eout | is.nan(tall)] = 0
  if (a > 0) {
    wp = which(shape > 1)
    tall[wp] = tall[wp] * ((4 * shape[wp])/(shape[wp] + 1)^2)^a
  }
  return(tall)
}

#' @export
#' @rdname scan_stat
stat.binom = function(yin, yout, ty, popin, popout, tpop) {
  py_in = popin - yin
  py_out = popout - yout
  
  tall = yin * (log(yin) - log(popin)) + 
    py_in * (log(py_in) - log(popin)) + 
    yout * (log(yout) - log(popout)) + 
    py_out * (log(py_out) - log(popout)) -
    ty * log(ty) - (tpop - ty) * log(tpop - ty) + 
    tpop * log(tpop)
  # correct test statistics for NaNs
  tall[yin/popin <= yout/popout | is.nan(tall)] = 0
  return(tall)
}

arg_check_scan_stat = function(yin, ty, ein, tpop, popin, type,
                               a, shape) {
  if (is.null(ein) & is.null(popin)) {
    stop("One of ein and popin must be provided")
  }
  if (!is.null(ein)) {
    if (length(yin) != length(ein)) {
      stop("yin and ein must be the same length")
    }
  }
  if (!is.null(popin)) {
    if (length(yin) != length(popin)) {
      stop("yin and popin must be the same length")
    }
  }
  if (length(type) != 1 | !is.element(type, c("poisson", "binomial"))) {
    stop("type must be poisson or binomial")
  }
  if (length(a) != 1 | a < 0 | a > 1) {
    stop("a must be a number between 0 and 1")
  }
  if (length(shape) != 1 & length(shape) != length(yin)) {
    stop("shape must be a single value or have length equal to length(yin)")
  }
  if (type == "binomial" & is.null(tpop)) {
    stop("tpop must be provided for the binomial model")
  }
}