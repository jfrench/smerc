#' Uses greedy algorithm to construct new candidate zones
#' 
#' \code{smst_zones_internal2} constructs a sequence of candidate zones from a starting region using a greedy algorithm that maximizes the scan statitistic subject to connectedness, population, and distance constraints.  It is an internal function meant to do the heavy lifting for the smst.zones function.  Each region is used as a starting region.
#' 
#' @param all_neighbors	A list with elements containing all neighbors for each region (in ascending order of distance from the the region).  A region is itself included in its neighbors.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param e The expected number of cases for each region.
#' @param ty The total number of cases.
#' @param max_pop The population upperbound for a candidate zone.
#' @param maxonly A logical value indicating whether only the maximum of the log likelihood ratio test statistic across all zones should be returned or more detailed information that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone.
#' @return Returns a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. Or the largest test statistic if maxonly = TRUE.
#' @author Joshua French
#' @importFrom matrixStats colMaxs


smst_max_zone()

smst_zones_internal2 = function(all_neighbors, cases, pop, w, e, ty, max_pop, maxonly = TRUE)
{
  stop = FALSE
  i = 1
  while(!stop)# && counter < 15)
  {
    # get current zone to extend
    c_zone = uz[[i]]
    
    # current max set of neighbors for current zone
    cmn = max_neighbors[[i]]
    
    # which potential regions are connected to c_zone
    # and satisfy the population constraints
    connected = cmn[which(matrixStats::colMaxs(w[c_zone, cmn, drop = FALSE]) == 1)]
    # new popin when adding each potential neighbor to c_zone
    p_popin = popin[i] + pop[connected]
    
    # candidate regions that satisfy constraints
    in_size = which(p_popin <= max_pop)
    cand_regions = connected[in_size]
    
    if(length(cand_regions) > 0)
    {
      # yin and ein for candidate zones (cur zone plus cand_regions)
      yin_cand = yin[i] + cases[cand_regions]
      ein_cand = ein[i] + e[cand_regions]
      
      # test statistics for candidate locations
      stat_cand = scan.stat(yin = yin_cand, ein = ein_cand,
                            eout = ty - ein_cand, ty = ty)
      
      # index of max stat_cand
      max_idx = which.max(stat_cand)
      
      # update for next best zone
      loglikrat[i + 1] = stat_cand[max_idx]
      yin[i + 1] = yin_cand[max_idx]
      ein[i + 1] = ein_cand[max_idx]
      popin[i + 1] = p_popin[in_size[max_idx]]
      # regions_used = c(regions_used, cand_regions[max_idx])
      max_neighbors[[i + 1]] = setdiff(max_neighbors[[i]], cand_regions[max_idx])
      uz[[i + 1]] = c(c_zone, cand_regions[max_idx])
      i = i + 1 # update counter
    }else # end algorithm
    {
      stop = TRUE
    }
  }
  # convert back to original location ids
  # return only the zones that meet constraint for population upper bound
  if(maxonly){
    return(max(loglikrat))
  }else
  {
    which_max = which.max(loglikrat)
    
    return(list(uz = uz[[which_max]],
                loglikrat = loglikrat[which_max],
                cases = yin[which_max],
                expected = ein[which_max],
                population = popin[which_max]))
  }
}


#' Determine sequence constrained greedy zones.
#' 
#' \code{smst_max_zone} finds the most likely zone using the maxima likelihood first scan algorithm.  It is an internal function meant to do the heavy lifting for the smst.zones function.
#' 
#' @param neighbors A vector containing the neighbors for a region (in ascending order of distance from the region).  The region itself is included in its neighbors.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param e The expected number of cases for each region.
#' @param ty The total number of cases.
#' @param max_pop The population upperbound for a candidate zone.
#' @param loglikrat is a vector to store the test statistic associated with each element of uz.  The first element already has the test statistic of uz[[1]]
#' @param yin is a vector to store the cases in each element of uz.  The first element already has the number of cases for uz[[1]]
#' @param ein is a vector to store the expected number of cases associated with each element of uz.  The first element already has the expected cases for uz[[1]]
#' @param popin s a vector to store the total population associated with each element of uz.  The first element already has the total population of uz[[1]]
#' @param type A character vector indicating whether only the maximum of the log likelihood ratio test statistic across all zones should be returned ("maxonly"), a more detailed list of information that includes the location id of each zone and the associated test statistic, counts, expected counts, and population for the zone with the max statistic  ("pruned"), or "all", which returns the same information for the sequence of candidate zones.
#' @return Returns If type = "all", a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. Or the largest test statistic if type = "maxonly".  Or the same information included in "all" for the max statistic zone if type = "pruned".
#' @author Joshua French
#' @importFrom matrixStats colMaxs
#' @examples
#' # load data
#' data(nydf)
#' data(nyw)
#' 
#' # create relevant data
#' coords = nydf[,c("longitude", "latitude")]
#' cases = floor(nydf$cases)
#' pop = nydf$population
#' w = nyw
#' ex = sum(cases)/sum(pop)*pop
#' ubpop = 0.5
#' ubd = 0.5
#' lonlat = TRUE
#' ty = sum(cases)   # total number of cases
#' # intercentroid distances
#' d = sp::spDists(as.matrix(coords), longlat = lonlat)
#' # upperbound for population in zone
#' max_pop = ubpop *sum(pop)
#' # upperbound for distance between centroids in zone
#' max_dist = ubd * max(d)
#' # create list of neighbors for each region (inclusive of region itself)
#' all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
#' # find the smst max zone
#' smst_max_zone(all_neighbors[[1]], cases, pop, w, ex, ty, max_pop)
#' smst_max_zone(all_neighbors[[1]], cases, pop, w, ex, ty, max_pop, "pruned")
#' bigout = smst_max_zone(all_neighbors[[1]], cases, pop, w, ex, ty, max_pop, "all")
#' head(bigout)

smst_max_zone = function(neighbors, cases, pop, w, e, ty, max_pop, type = "maxonly")
{
  
  loglikrat = yin = ein = popin = numeric(length(neighbors))
  region = neighbors[1]
  yin[1] = cases[region]
  ein[1] = e[region]
  popin[1] = pop[region]
  loglikrat[1] = scan.stat(yin[1], ein[1], ty - ein[1], ty)
  uz = max_neighbors = vector("list", length(neighbors))
  uz[[1]] = region
  max_neighbors[[1]] = neighbors[-1]

  stop = FALSE
  i = 1
  while(!stop)# && counter < 15)
  {
    # get current zone to extend
    c_zone = uz[[i]]
    
    # current max set of neighbors for current zone
    cmn = max_neighbors[[i]]
    
    # which potential regions are connected to c_zone
    # and satisfy the population constraints
    connected = cmn[which(matrixStats::colMaxs(w[c_zone, cmn, drop = FALSE]) == 1)]
    # new popin when adding each potential neighbor to c_zone
    p_popin = popin[i] + pop[connected]
    
    # candidate regions that satisfy constraints
    in_size = which(p_popin <= max_pop)
    cand_regions = connected[in_size]
    
    if(length(cand_regions) > 0)
    {
      # yin and ein for candidate zones (cur zone plus cand_regions)
      yin_cand = yin[i] + cases[cand_regions]
      ein_cand = ein[i] + e[cand_regions]
      
      # test statistics for candidate locations
      stat_cand = scan.stat(yin = yin_cand, ein = ein_cand,
                            eout = ty - ein_cand, ty = ty)
      
      # index of max stat_cand
      max_idx = which.max(stat_cand)
      
      # update for next best zone
      loglikrat[i + 1] = stat_cand[max_idx]
      yin[i + 1] = yin_cand[max_idx]
      ein[i + 1] = ein_cand[max_idx]
      popin[i + 1] = p_popin[in_size[max_idx]]
      # regions_used = c(regions_used, cand_regions[max_idx])
      max_neighbors[[i + 1]] = setdiff(max_neighbors[[i]], cand_regions[max_idx])
      uz[[i + 1]] = c(c_zone, cand_regions[max_idx])
      i = i + 1 # update counter
    }else # end algorithm
    {
      stop = TRUE
    }
  }
  # return desired output
  if(type == "maxonly"){
    return(max(loglikrat))
  }else if (type == "pruned")
  {
    which_max = which.max(loglikrat)
    
    return(list(locids = uz[[which_max]],
                loglikrat = loglikrat[which_max],
                cases = yin[which_max],
                expected = ein[which_max],
                population = popin[which_max]))
  }else if(type == "all")
  {
    return(list(locids = uz,
                loglikrat = loglikrat,
                cases = yin,
                expected = ein,
                population = popin))
  }
}

