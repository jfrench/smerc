#' Determine zones for flexibly shaped spatial scan test
#' 
#' \code{flex.zones} determines the unique zones to consider for the flexibly shaped spatial scan test of Tango and Takahashi (2005).  The algorithm uses a breadth-first search to find all subgraphs connected to each vertex (region) in the data set of size \eqn{k} or less.  
#' 
#' @param coords The centroid coordinates for each region.
#' @param w The binary spatial adjacency matrix of the regions.
#' @param k The maximum number of regions to include in a zone.
#' @param lonlat A logical indioating whether the coordinates are longitude/latitude.  If so, the great circle distance is used in computing the nearest/neighbor distance matrix.
#' @param parallel A logical value indicating whether the \code{parallel} package should be used to speedup computations.
#' @return Returns a list of zones to consider for clustering.  Each element of the list contains a vector with the location ids of the regions in that zone.
#' @author Joshua French
#' @importFrom spdep knearneigh
#' @importFrom parallel mclapply
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' flex.zones(coords = coords, w = nyw, k = 3, lonlat = TRUE)
#' 
flex.zones = function(coords, w, k = 10, lonlat = FALSE, parallel = TRUE)
{
  N = nrow(coords)
  
  mynn = cbind(1:N, spdep::knearneigh(as.matrix(coords), k = (k - 1), longlat = lonlat)$nn)

  fcall = lapply
  if(parallel) fcall = parallel::mclapply
  fcall_list = list(X = as.list(1:N), function(i)
  {
    connected_subgraphs(w = w[mynn[i, ], mynn[i, ]],
                        nn = mynn[i, ], k = k)
  })
  
  czones = unlist(do.call(fcall, fcall_list), 
                  use.names = FALSE, 
                  recursive = FALSE)
  
  if(parallel)
  {
    return(unique(parallel::mclapply(czones, sort)))
  }else
  {
    return(unique(lapply(czones, sort)))
  }
}

# takes a spatial adjacency matrix and the 
# index of the locations in the spatial adjacency
# matrix
connected_subgraphs = function(w, nn, k)
{
  # storage list, of length k
  listi = vector("list", k)
  listi[[1]] = as.list(1)
  
  # index of neighbors for each region
  nbi = apply(w, 2, function(x) which(x!=0))
  
  for(j in 2:k)
  {
    if(!is.null(listi[[j - 1]]))
    {
     newset = unique(unlist(
      lapply(listi[[j - 1]], function(r)
      {
        sd = setdiff(unique(unlist(nbi[r])), r)
        if(length(sd) > 0)
        sapply(sd,
               function(u) sort(c(r, u)), 
               simplify = FALSE)
      }), recursive = FALSE))
     if(length(newset) > 0)
     {
       listi[[j]] = newset
     }else
     {
       j = k + 1
     }
    }
  }  
  
  return(sapply(unlist(listi, recursive = FALSE), function(l) nn[l]))
}