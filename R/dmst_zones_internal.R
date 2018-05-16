## Determine zones using the dynamic minimum spanning tree scan test of Assuncao et al. (2006) 
## 
## \code{dmst_zones_internal} determines the zones that produce the largest test statistic using a greedy algorithm.  Specifically, starting individually with each region as a starting zone, new (connected) regions are added to the current zone in the order that results in the largest test statistic.  This is used to implement the dynamic minimum spanning tree (dmst) scan test of Assuncao et al. (2006).  See also the maxima likelihood first scan method of Yao et al. (2011).
## 
## Every zone considered must have a total population less than \code{ubpop * the total population across all regions} in the study area.  Additionally, the maximum intercentroid distance for the regions within a zone must be no more than \code{ubd * the maximum intercentroid distance across all regions}. 
## 
## @param all_neighbors	A list of length \eqn{n} containing the neighbors of each region (inclusive) subject to the distance constraint.
## @param cases The number of cases in each region.
## @param pop The population size of each region.
## @param w The binary spatial adjacency matrix.
## @param ex The expected number of cases for each region. The default is calculated under the constant risk hypothesis.
## @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
## @param parallel A logical value indicating whether finding candidate zones should use multiple cores.  The default is \code{FALSE}.  Multiple cores are used if \code{paralle = TRUE}.
## @param maxonly A logical value indicating whether to return only the maximum test statistic across all candidate zones.  Default is \code{FALSE}.

## @return Returns a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. If \code{maxonly = TRUE}, then only the maximum test statistic across all of these zones is returned.
## @author Joshua French
## @importFrom utils tail
## @references Assuncao, R.M., Costa, M.A., Tavares, A. and Neto, S.J.F. (2006). Fast detection of arbitrarily shaped disease clusters, Statistics in Medicine, 25, 723-742.
## 
## Yao, Z., Tang, J., & Zhan, F. B. (2011). Detection of arbitrarily-shaped clusters using a neighbor-expanding approach: A case study on murine typhus in South Texas. International Journal of Health Geographics, 10(1), 1.
## @export
## @examples 
## data(nydf)
## data(nyw)
## coords = as.matrix(nydf[,c("longitude", "latitude")])
## # obtain distance matrix of intercentroid distances
## d = sp::spDists(as.matrix(coords), longlat = TRUE)
## # upperbound for distance between centroids in zone
## max_dist = 0.5 * max(d)
## # find all neighbors from each starting zone within distance upperbound
## all_neighbors = lapply(seq_along(nydf$cases), function(i) which(d[i,] <= max_dist))
## # expected counts
## ex = sum(nydf$cases)/sum(nydf$pop)*nydf$pop
## # find zone with max statistic starting from each individual region
## max_zones = dmst_zones_internal(all_neighbors, cases = floor(nydf$cases), 
##                                 nydf$pop, w = nyw, ubpop = 0.5)
## head(max_zones)
## # return only the largest test statistic across all the candidate zones
## max_zones2 = dmst.zones(coords, cases = floor(nydf$cases), 
##                         nydf$pop, w = nyw, ubpop = 0.5, 
##                         ubd = .5, lonlat = TRUE, maxonly = TRUE)

# data(nydf)
# data(nyw)
# coords = nydf[,c("longitude", "latitude")]
# cases = floor(nydf$cases)
# pop = nydf$population
# w = nyw
# e = sum(cases)/sum(pop)*pop
# ubpop = 0.5
# ubd = 0.5
# lonlat = TRUE
# parallel = FALSE 
# maxonly = TRUE
dmst_zones_internal = function(all_neighbors, cases, pop, w, ex, ubpop = 0.5, cl = NULL, maxonly = FALSE)
{
  # sanity checking
  arg_check_dmst_zones_internal(all_neighbors, cases, pop, w, ex, ubpop, FALSE)
  
  # setup various arguments and such
  ty = sum(cases)   # total number of cases
  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # should only max be returned, or a pruned version
  type = ifelse(maxonly, "maxonly", "pruned")
  
  # set up finding zones with max stat from each starting region
  fcall = pbapply::pblapply
  fcall_list = list(X = as.list(seq_along(all_neighbors)), FUN = function(i){
    # find zone with max stat, starting from each region
    dmst_max_zone(i, all_neighbors[[i]], cases, pop, w, ex, ty, max_pop, type = type)
  }, cl = cl)

  # obtain list of zones with maximum statistic from each starting region (or the max statistic)
  max_zones = do.call(fcall, fcall_list)
  
  # return max statistic or the max statistic zones
  if(maxonly){
    return(max(unlist(max_zones)))
  }else{
    return(max_zones)
  }
}

arg_check_dmst_zones_internal = function(all_neighbors, cases, pop, w, ex, ubpop, parallel)
{
  if(!is.list(all_neighbors)) stop("all_neighbors should be a list")
  if(length(all_neighbors) != length(cases)) stop("length(all_neighbors) != length(cases)")
  if(length(cases) != length(pop)) stop('length(cases) != length(pop)')
  if(length(cases) != nrow(w)) stop('length(cases) != nrow(w)')
  if(length(cases) != length(ex)) stop('length(cases) != length(ex)')
  if(length(ubpop) != 1 | !is.numeric(ubpop)) stop("ubpop should be a single number")
  if(ubpop <= 0 | ubpop > 1) stop("ubpop not in (0, 1]")
  if(length(parallel) != 1 || !is.logical(parallel)) stop("parallel must be a single logical value")
}