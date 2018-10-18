#' Spatial scan statistic
#' 
#' \code{scan.stat} calculates the spatial scan statistic 
#' for a zone (set of spatial regions).
#' 
#' @param yin The total number of cases in the zone.
#' @param ein The expected number of cases in the zone.  
#' Conventionally, this is the estimated overall 
#'   disease risk across the study area, multiplied by the 
#'   total population size of the zone.
#' @param eout The expected number of cases outside the 
#'   zone.  This should be \code{ty - ein}, but
#'   this is not done for efficiency purposes.
#' @param ty The total number of cases in the study area.
#' @param type The type of scan statistic to implement. 
#'   Currently, only \code{"poisson"} is implemented.
#' @param shape The shape of the ellipse, which is the ratio
#' of the length of the longest and shortest axes of the ellipse.  The
#' default is 1, meaning it is a circle.
#' @param a A tunning parameter for the adjusted log-likelihood ratio.  See details.

#' @return A vector of scan statistics.
#' @author Joshua French
#' @export
#' @references 
#' Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' # statistic for most likely cluster of 
#' # New York leukemia data
#' scan.stat(106, 62.13, 552 - 62.13, 552)
scan.stat = function(yin, ein, eout, ty, type = "poisson", a = 0, shape = 1) {
  if (length(yin) != length(ein)) {
    stop("length(yin) != length(ein)")
  }
  if (length(yin) != length(eout)) {
    stop("length(yin) != length(eout)")
  }
  if (type == "poisson") {
    yout = ty - yin
    tall = yin * (log(yin) - log(ein)) + yout * (log(yout) - log(eout))
    # correct test statistics for NaNs
    tall[yin/ein <= yout/eout | is.nan(tall)] = 0
    if (a > 0) {
      wp = which(shape > 1)
      tall[wp] = tall[wp] * ((4 * shape[wp])/(shape[wp] + 1)^2)^a
    }
  }
  return(tall)
}