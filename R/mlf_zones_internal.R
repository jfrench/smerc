#' Determine sequence constrained greedy zones.
#' 
#' \code{mlf_zones_internal} finds the most likely zone using the maxima likelihood first scan algorithm.  It is an internal function meant to do the heavy lifting for the smst.zones function.
#' 
#' @param uz	A list of unique zones from the algorithm.  The first element is the most likely cluster for a single region.
#' @param max_neighbors A list.  The first element has the set of potential regions to add to the zone for uz[[1]].
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
#' @param maxonly A logical value indicating whether only the maximum of the log likelihood ratio test statistic across all zones should be returned or more detailed information that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone.
#' @return Returns a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. Or the largest test statistic if maxonly = TRUE.
#' @author Joshua French
#' @importFrom matrixStats colMaxs

mlf_zones_internal = function(uz, max_neighbors, cases, pop, w, e, ty, max_pop, loglikrat, yin, ein, popin, maxonly = TRUE)
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

# cog_zones_internal = function(all_neighbors, cases, pop, w, e, max_pop, ty, maxonly = TRUE)
# {
#   N = length(cases) # number of regions
#   
#   # scan statistics for each individual region
#   stat = scan.stat(cases, e, ty - e, ty) * (pop <= max_pop)
#   
#   # regions that can't be a zone or will never increase the likelihood of a zone are "used", i.e., can't be considered for inclusion in a zone
#   regions_used = which(pop > max_pop | cases == 0)
#   regions_left = setdiff(1:N, regions_used)
#   
#   # for storage
#   loglikrat = yin = ein = popin = numeric(length(regions_left))
#   
#   # vector of unique zones, potential neighbors for unique zones based on distance
#   uz =  max_neighbors = vector("list", length(regions_left))
#   
#   # find starting zone (the region with the largest stat)
#   which_max_stat = which.max(stat[regions_left])
#   uz[[1]] = regions_left[which_max_stat]
#   
#   # determine the regions that have been used in a zone
#   regions_used = c(uz[[1]], regions_used)
#   # determine the regions that have not been used in a zone
#   regions_left = setdiff(1:N, regions_used)
#   
#   # current zones that are being added to
#   cz = 1
#   # update test statistic, number of cases, expected number
#   # of cases, total population for zone
#   loglikrat[1] = stat[uz[[1]]]
#   yin[1] = cases[uz[[1]]]
#   ein[1] = e[uz[[1]]]
#   popin[1] = pop[uz[[1]]]
#   max_neighbors[[1]] = setdiff(all_neighbors[[uz[[1]]]], regions_used)
#   
#   counter = 2
#   while(length(regions_used) < N)# && counter < 15)
#   {
#     cur_best = -1 # current best statistic
#     idx_best = NULL # index of the best statistic
#     zone_best = NULL # the best zone of the choice considered
#     cmn_best = NULL # the remaining common neighbors to consider for the best zone
#     yin_best = NULL
#     ein_best = NULL
#     popin_best = NULL
#     
#     cz_active = !logical(length(cz)) # all current zones are active
#     which_cz_best = 0 # index of cz that produces maximum statistic when adding adjacent region
#     for(i in seq_along(cz)) # for each active zone
#     {
#       # locids of the current active zone
#       c_zone = uz[[cz[i]]]
#       
#       # current max set of neighbors for current active zone
#       cmn = max_neighbors[[cz[[i]]]]
#       
#       # which potential regions are connected to c_zone
#       # and satisfy the population constraints
#       connected = cmn[which(matrixStats::colMaxs(w[uz[[cz[i]]], cmn, drop = FALSE]) == 1)]
#       # new popin when adding each potential neighbor to c_zone
#       p_popin = popin[[cz[i]]] + pop[connected]
#       
#       # candidate regions that satisfy constraints
#       in_size = which(p_popin <= max_pop)
#       cand_regions = connected[in_size]
#       
#       stat_cand = -2
#       if(length(cand_regions) > 1)
#       {
#         # yin and ein for candidate zones (cur zone plus cand_regions)
#         yin_cand = yin[cz[i]] + cases[cand_regions]
#         ein_cand = ein[cz[i]] + e[cand_regions]
#         
#         # test statistics for candidate locations
#         stat_cand = scan.stat(yin = yin_cand, ein = ein_cand,
#                               eout = ty - ein_cand, ty = ty)
#         
#         # index of max stat_cand
#         max_idx = which.max(stat_cand)
#         
#         # update for best zone
#         stat_cand = stat_cand[max_idx]
#         zone_cand = c(uz[[cz[i]]], cand_regions[max_idx])
#         yin_cand = yin_cand[max_idx]
#         ein_cand = ein_cand[max_idx]
#         popin_cand = p_popin[in_size[max_idx]]
#         cmn_cand = cand_regions[-max_idx] 
#       }
#       
#       if(stat_cand > cur_best)
#       {
#         cur_best = stat_cand
#         zone_best = zone_cand
#         yin_best = yin_cand
#         ein_best = ein_cand
#         popin_best = popin_cand
#         which_cz_best = i
#         cmn_best = cmn_cand
#       }
#       # if no connected regions can be added to 
#       # current zone, deactivate the zone
#       if(stat_cand == -2) cz_active[i] = FALSE
#     }
#     
#     # determine if any new zones should be started
#     # i.e., the statistic for a single region
#     # is more than the statistic you get from combining
#     # a potential region with the active zones
#     max_stat_left = max(stat[regions_left])
#     if(max_stat_left > cur_best)
#     {
#       # determine which index in stat[regions_left] has the maximum
#       wmsl = which.max(stat[regions_left])
#       zone_best = regions_left[wmsl]
#       cur_best = max_stat_left
#       yin_best = cases[zone_best]
#       ein_best = e[zone_best]
#       popin_best = pop[zone_best]
#       # include potential neighbors, 
#       # minus the best region itself and regions left
#       cmn_best = setdiff(all_neighbors[[zone_best]], c(zone_best, regions_left))
#     }
#     
#     uz[[counter]] = zone_best
#     max_neighbors[[counter]] = cmn_best
#     loglikrat[counter] = cur_best
#     yin[counter] = yin_best
#     ein[counter] = ein_best
#     popin[counter] = popin_best
#     
#     if(length(zone_best) == 1)
#     {
#       # add new zone to list
#       cz = c(cz, counter)
#     }else
#     {
#       # add new zone to list
#       cz_best = cz[which_cz_best]
#       cz_delete = cz[!cz_active]
#       cz = setdiff(c(cz, counter), c(cz_delete, cz_best))
#     }
#     regions_used = c(regions_used, tail(zone_best, 1))
#     regions_left = setdiff(regions_left, tail(zone_best, 1))
#     counter = counter + 1
#   }
#   # convert back to original location ids
#   # return only the zones that meet constraint for population upper bound
#   if(maxonly){
#     return(max(loglikrat))
#   }else
#   {
#     return(list(zones = uz,
#                 loglikrat = loglikrat,
#                 cases = yin,
#                 expected = ein,
#                 population = popin))
#   }
# }