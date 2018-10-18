#' Determine zones for flexibly shaped spatial scan test
#'
#' \code{rflex.zones} determines the unique zones to
#' consider for the flexibly shaped spatial scan test of
#' Tango and Takahashi (2012).  The algorithm uses a
#' breadth-first search to find all subgraphs connected to
#' each vertex (region) in the data set of size \eqn{k} or
#' less with the constraint that the middle p-value of each
#' region must be less than \code{alpha1}.
#'
#' @inheritParams flex.test
#' @param nn An n by k matrix providing the k nearest
#'   neighbors of each region, presumably produced by the
#'   \code{\link{knn}} function.
#' @param alpha1 The middle p-value threshold.
#' @param progress A logical value indicating whether the
#'   progress of constructing the zones should be reported.
#'   The default is FALSE.
#' @return Returns a list of zones to consider for
#'   clustering.  Each element of the list contains a vector
#'   with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Tango, T. and Takahashi, K. (2012), A
#'   flexible spatial scan statistic with a restricted
#'   likelihood ratio for detecting disease clusters.
#'   Statist. Med., 31: 4207-4218. <doi:10.1002/sim.5478>
#' @seealso rflex.midp
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$longitude, nydf$latitude)
#' nn = knn(coords, longlat = FALSE, k = 50)
#' cases = floor(nydf$cases)
#' ex = nydf$pop * sum(cases)/sum(nydf$pop)
#' zones = rflex.zones(nn, w = nyw, cases = cases, ex = ex)
rflex.zones = function(nn, w, cases, ex, 
                       alpha1 = 0.2,
                       cl = NULL, 
                       progress = FALSE) {
  N = nrow(nn)
  
  # compute mid p-value 
  p = rflex.midp(cases, ex)
  keep = which(p < alpha1)
  remove = setdiff(seq_len(N), keep)
  
  # remove connections when p >= alpha1  
  w[,remove] = 0
  
  if (progress) {
    message("constructing connected subgraphs:")
    fcall = pbapply::pblapply
  } else {
    fcall = lapply
  }

  # allcg = vector("list", length(keep))
  # j = 1
  # for (i in keep) {
  #   idxi = intersect(nn[i, ], c(which(w[i, ] == 1), i))
  #   allcg[[j]] = scsg(idxi, w[,idxi, drop = FALSE])
  #   j = j + 1
  # }
  fcall_list = list(X = keep, function(i, ...) {
    idxi = intersect(nn[i, ], keep)
    scsg(idxi, w[,idxi, drop = FALSE])
  }, cl = cl)  
  czones = unlist(do.call(fcall, fcall_list), 
                  use.names = FALSE, 
                  recursive = FALSE)
  czones[distinct(czones)]
}

