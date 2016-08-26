## Determine sequence of constrained greedy zones.
## 
## \code{dmst_max_zone} finds the sequence of connected regions that maximize the likelihood from a starting region.  The sequence continues until no region can be added to the current zone due to constraints.
## @param start The region to start extending the zone from.
## @param neighbors A vector containing the neighbors for a region (in ascending order of distance from the region).  The region itself is included in its neighbors.
## @param cases The number of cases in each region.
## @param pop The population size of each region.
## @param w The binary spatial adjacency matrix.
## @param ex The expected number of cases for each region.
## @param ty The total number of cases.
## @param max_pop The population upperbound for a candidate zone.
## @param type A character vector indicating whether only the maximum of the log likelihood ratio test statistic across all zones should be returned ("maxonly"), a more detailed list of information that includes the location id of each zone and the associated test statistic, counts, expected counts, and population for the zone with the max statistic  ("pruned"), or "all", which returns the same information for the sequence of candidate zones.
## @return Returns If type = "all", a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. Or the largest test statistic if type = "maxonly".  Or the same information included in "all" for the max statistic zone if type = "pruned".
## @author Joshua French
## @examples
## # load data
## data(nydf)
## data(nyw)
## 
## # create relevant data
## coords = nydf[,c("longitude", "latitude")]
## cases = floor(nydf$cases)
## pop = nydf$population
## w = nyw
## ex = sum(cases)/sum(pop)*pop
## ubpop = 0.5
## ubd = 0.5
## lonlat = TRUE
## ty = sum(cases)   # total number of cases
## # intercentroid distances
## d = sp::spDists(as.matrix(coords), longlat = lonlat)
## # upperbound for population in zone
## max_pop = ubpop *sum(pop)
## # upperbound for distance between centroids in zone
## max_dist = ubd * max(d)
## # create list of neighbors for each region (inclusive of region itself)
## all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
## # find the dmst max zone
## dmst_max_zone(start = 1, all_neighbors[[1]], cases, pop, w, ex, ty, max_pop)
## dmst_max_zone(start = 1, all_neighbors[[1]], cases, pop, w, ex, ty, max_pop, "pruned")
## dmst_max_zone(start = 2, all_neighbors[[2]], cases, pop, w, ex, ty, max_pop, "pruned")
## bigout = dmst_max_zone(start = 1, all_neighbors[[1]], cases, pop, w, ex, ty, max_pop, "all")
## head(bigout)
dmst_max_zone = function(start, neighbors, cases, pop, w, ex, ty, max_pop, type = "maxonly")
{
  loglikrat = yin = ein = popin = numeric(length(neighbors))
  region = start
  yin[1] = cases[region]
  ein[1] = ex[region]
  popin[1] = pop[region]
  loglikrat[1] = scan.stat(yin[1], ein[1], ty - ein[1], ty)
  uz = max_neighbors = vector("list", length(neighbors))
  uz[[1]] = region
  max_neighbors[[1]] = setdiff(neighbors, region)

  stop = FALSE
  # double check that pop constraint is satisfied
  # if not, stop and set loglikrat to 0 (so that it's not the maximum in the earch)
  if(popin[1] > max_pop){
    stop = TRUE
    loglikrat[1] = 0
  }
  
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
      ein_cand = ein[i] + ex[cand_regions]
      
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

