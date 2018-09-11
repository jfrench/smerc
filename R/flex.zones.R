#' Determine zones for flexibly shaped spatial scan test
#' 
#' \code{flex.zones} determines the unique zones to consider for the flexibly shaped spatial scan test of Tango and Takahashi (2005).  The algorithm uses a breadth-first search to find all subgraphs connected to each vertex (region) in the data set of size \eqn{k} or less.  
#' 
#' @inheritParams flex.test
#' @param progress A logical value indicating whether a progress bar should be displayed.  The default is \code{TRUE}.
#' @return Returns a list of zones to consider for clustering.  Each element of the list contains a vector with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' zones = flex.zones(coords, w = nyw, k = 3, longlat = TRUE)
flex.zones = function(coords, w, k = 10, longlat = FALSE, 
                      cl = NULL, progress = TRUE) {
  N = nrow(coords)
  
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  mynn = t(apply(d, 1, order)[seq_len(k), ])

  if (progress) {
    message("constructing connected subgraphs:")
    fcall = pbapply::pblapply
  } else {
    fcall = lapply
  }

  fcall_list = list(X = as.list(1:N), function(i, ...) {
    connected_subgraphs(w = w[mynn[i, ], mynn[i, ]],
                        nn = mynn[i, ], k = k)
  }, cl = cl)  
  czones = unlist(do.call(fcall, fcall_list), 
                  use.names = FALSE, 
                  recursive = FALSE)
  czones[distinct(czones)]
}

# takes a spatial adjacency matrix and the 
# index of the locations in the spatial adjacency
# matrix
connected_subgraphs = function(w, nn, k) {
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
