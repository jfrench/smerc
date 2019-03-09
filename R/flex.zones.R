#' Determine zones for flexibly shaped spatial scan test
#' 
#' \code{flex.zones} determines the unique zones to consider for the flexibly shaped spatial scan test of Tango and Takahashi (2005).  The algorithm uses a breadth-first search to find all subgraphs connected to each vertex (region) in the data set of size \eqn{k} or less.  
#' 
#' @inheritParams flex.test
#' @inheritParams rflex.zones
#' @param progress A logical value indicating whether a progress bar should be displayed.  The default is \code{TRUE}.
#' @return Returns a list of zones to consider for clustering.  Each element of the list contains a vector with the location ids of the regions in that zone.
#' @author Joshua French
#' @export
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = cbind(nydf$x, nydf$y)
#' zones = flex.zones(coords, w = nyw, k = 3)
#' \dontrun{
#' # see what happens when verbose = TRUE
#' zones = flex.zones(coords, w = nyw, k = 3, verbose = TRUE)
#' }
flex.zones = function(coords, w, k = 10, longlat = FALSE,
                      cl = NULL, progress = TRUE, verbose = FALSE) {
  nn = knn(coords = coords, longlat = longlat, k = k)
  N = nrow(coords)

  if (!verbose) {
    fcall_list = list(X = seq_len(N), function(i, ...) {
      scsg(nn[[i]], w[, nn[[i]], drop = FALSE])
    }, cl = cl)

    # determine which apply function to use
    if (progress) {
      message("constructing connected subgraphs:")
      fcall = pbapply::pblapply
    } else {
      fcall = lapply
    }

    # determine zones
    czones = unlist(do.call(fcall, fcall_list),
                    use.names = FALSE, recursive = FALSE)
  } else {
    czones = list()
    for (i in seq_len(N)) {
      if (verbose) {
        message(paste(i, "/", N, ". Starting region ", i,
                      " at ", Sys.time(), ".", sep = ""))
      }
      czones = combine.zones(czones,
                             scsg(nn[[i]], w[, nn[[i]], drop = FALSE]))
    }
    return(czones)
  }
  # determine distinct zones
  czones[distinct(czones)]
}