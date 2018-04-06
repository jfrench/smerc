#' Determine the candidate zone using the maxima likelihood first algorithm of Yao et al. (2011).
#' 
#' \code{mlf.zones} determines the most likely cluster zone obtained by implementing the maxima likelihood first scann method of Yao et al. (2011).  Note that this is really just a special case of the dynamic minimum spanning tree (SMST) algorithm of Assuncao et al. (2006) 
#' 
#' Each step of the mlf scan test seeks to maximize the likelihood ratio test statistic used in the original spatial scan test (Kulldorff 1997).  The first zone considered is the region that maximizes this likelihood ration test statistic, providing that no more than \code{ubpop} proportion of the total population is in the zone.  The second zone is the first zone and the connected region that maximizes the scan statistic, subject to the population and distance constraints.  This pattern continues until no additional zones can be added due to population or distance constraints.
#' 
#' Every zone considered must have a total population less than \code{ubpop * sum(pop)} in the study area.  Additionally, the maximum intercentroid distance for the regions within a zone must be no more than \code{ubd * the maximum intercentroid distance across all regions}. 
#' 
#' @inheritParams dmst.zones
#' @param type One of \code{"maxonly"}, \code{"pruned"}, or \code{"all"}. Specifying \code{"maxonly"} returns only the maximum test statistic across all candidate zones, \code{"pruned"} returns information for the zone with the largest test statistic, while \code{"all"} returns information for all candidate zones.  Default is \code{"pruned"}.
#' @return Returns a list that includes the location id of the zone and the associated test statistic, counts, expected counts, and population in the zone.  If \code{type = "all"}, then each of these elements is a list or vector corresponding to each respective candidate zone.
#' @author Joshua French
#' @references Yao, Z., Tang, J., & Zhan, F. B. (2011). Detection of arbitrarily-shaped clusters using a neighbor-expanding approach: A case study on murine typhus in South Texas. International journal of health geographics, 10(1), 1.
#' 
#' Assuncao, R.M., Costa, M.A., Tavares, A. and Neto, S.J.F. (2006). Fast detection of arbitrarily shaped disease clusters, Statistics in Medicine, 25, 723-742.
#' @export
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("x", "y")])
#' mlf.zones(coords, cases = floor(nydf$cases), pop = nydf$pop, w = nyw, lonlat = TRUE)

# data(nydf)
# data(nyw)
# coords = nydf[,c("longitude", "latitude")]
# cases = floor(nydf$cases)
# pop = nydf$population
# w = nyw
# ex = sum(cases)/sum(pop)*pop
# ubpop = 0.5
# ubd = 0.5
# lonlat = TRUE
# parallel = FALSE 
# maxonly = TRUE
# mlf.zones(coords, cases, pop, w, ex, lonlat = TRUE)
# mlf.zones(coords, cases, pop, w, ex, lonlat = TRUE, type = "maxonly")
# mlf.zones(coords, cases, pop, w, ex, lonlat = TRUE, type = "all")

mlf.zones = function(coords, cases, pop, w, ex = sum(cases)/sum(pop)*pop, ubpop = 0.5, ubd = 1, lonlat = FALSE, parallel = TRUE, type = "pruned")
{
  # sanity checking
  arg_check_mlf_zones(coords, cases, pop, w, ex, ubpop, ubd, lonlat, parallel, type)
  
  # total number of cases
  ty = sum(cases)
  
  # calculate test statistics for each individual region
  # yin, ein, eout, ty
  tobs = scan.stat(cases, ex, ty - ex, ty)
  
  # determine starting region for maxima likelihood first algorithm
  start = which.max(tobs)
  
  # intercentroid distances 
  d = sp::spDists(as.matrix(coords), longlat = lonlat)

  # upperbound for population in zone
  max_pop = ubpop *sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)

  # find neighbors of starting region
  start_neighbors = which(d[start,] <= max_dist)
  
  # return sequence of candidate zones (or a subset depending on type)
  dmst_max_zone(start, start_neighbors, cases, pop, w, ex, ty, max_pop, type)
}

arg_check_mlf_zones = function(coords, cases, pop, w, ex, ubpop, ubd, lonlat, parallel, type)
{
  if(ncol(coords) != 2) stop("coords should have 2 columns")
  if(nrow(coords) != length(cases)) stop("nrow(coords) != length(cases)")
  if(length(cases) != length(pop)) stop('length(cases) != length(pop)')
  if(length(cases) != nrow(w)) stop('length(cases) != nrow(w)')
  if(nrow(w) != ncol(w)) stop("w must be a square matrix")
  if(length(cases) != length(ex)) stop('length(cases) != length(ex)')
  if(length(ubpop) != 1 | !is.numeric(ubpop)) stop("ubpop should be a single number")
  if(ubpop <= 0 | ubpop > 1) stop("ubpop not in (0, 1]")
  if(length(ubd) != 1 | !is.numeric(ubd)) stop("ubd should be a single number")
  if(ubd <= 0 | ubd > 1) stop("ubd not in (0, 1]")
  if(length(lonlat) != 1 || !is.logical(lonlat)) stop("lonlat must be a single logical value")
  if(length(parallel) != 1 || !is.logical(parallel)) stop("parallel must be a single logical value")
  if(!is.element(type, c("maxonly", "pruned", "all"))) stop("invalid type")
}