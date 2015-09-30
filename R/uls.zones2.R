#' Determine sequence of ULS zones.
#' 
#' \code{uls.zones} determines the unique zones obtained by implementing the ULS (Upper Level Set) method of Patil and Taillie (2004).
#' 
#' The zones returned must have a total population less than ubpop * the total population of all regions in the study area.
#' 
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
#' @return Returns a list of zones to consider for clustering.  Each element of the list contains a vector with the location ids of the regions in that zone.
#' @author Joshua French
#' @import igraph
#' @references Patil, G. P., and Taillie, C. (2004). Upper level set scan statistic for detecting arbitrarily shaped hotspots. Environmental and Ecological Statistics, 11(2), 183-197.
#' @examples 
#' data(nydf)
#' data(nyw)
#' uls.zones2(cases = nydf$cases, pop = nydf$population, w = nyw)
uls.zones2 = function(cases, pop, w, ubpop = 0.5)
{
  # order rates from largest to smallest
  or = order(cases/pop, decreasing = TRUE)
  # reorder rows and columns by order of r
  w = w[or, ]
  w = w[, or]
  zones_list = vector("list", nrow(w))
  current_zones = NULL
  g = igraph::make_empty_graph(n = 0)
  for(i in 1:nrow(w))
  {
    # add new vertex
    g = igraph::add_vertices(g, nv = 1)
    # determine if there are any potential edges
    pot_edge = which(w[i, 1:i, drop = FALSE] == 1)
    # if there are, add them
    if(length(pot_edge) > 0) g = igraph::add_edges(g, c(rbind(i, pot_edge)))
    # c_list[[i]]$g = g
    # c_list[[i]]$cl = clusters(g)
    # zones_list[[i]] = tapply(1:i, c_list[[i]]$cl$membership, function(x) return(x))
    # return all zones for this threshold
    zones_list[[i]] = tapply(1:i, igraph::clusters(g)$membership, function(x) return(x))
  }
  # unique zones
  uz = unique(unlist(zones_list, recursive = FALSE))
  # convert back to original ids
  uz = lapply(uz, function(x) or[x])
  # upper bound for population in window
  ub = sum(pop) * ubpop
  # return only the zones that are unique
  popin = lapply(uz, function(x) sum(pop[x]))
  return(uz[which(popin <= ub)])
}