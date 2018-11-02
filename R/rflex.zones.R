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
#' @inheritParams rflex.test
#' @param nn An n by k matrix providing the k nearest
#'   neighbors of each region, presumably produced by the
#'   \code{\link{knn}} function.
#' @param pop The population size associated with each
#'   region.  The default is \code{NULL} since this argument
#'   is only needed for \code{type = "binomial"}.
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
#' pop = nydf$pop
#' ex = pop * sum(cases)/sum(pop)
#' # zones for poisson model
#' pzones = rflex.zones(nn, w = nyw, cases = cases, ex = ex)
#' # zones for binomial model
#' bzones = rflex.zones(nn, w = nyw, cases = cases, ex = ex,
#'                      type = "binomial", pop = pop)
rflex.zones = function(nn, w, cases, ex, alpha1 = 0.2, 
                       type = "poisson", pop = NULL, 
                       cl = NULL, progress = FALSE) {
  arg_check_rflex_zones(nn, w, cases, ex, alpha1, type, pop, 
                        progress)
 
  # compute mid p-value 
  p = rflex.midp(cases, ex, type = type, pop = pop)
  # determine which regions are "hot" (keep) or "cold" (remove)
  keep = which(p < alpha1)
  remove = setdiff(seq_along(ex), keep)
  
  # remove connections when p >= alpha1  
  w[,remove] = 0
  
  fcall_list = list(X = keep, function(i, ...) {
    idxi = intersect(nn[i, ], keep)
    scsg(idxi, w[,idxi, drop = FALSE])
  }, cl = cl)  
  
  # determine which apply function to use
  if (progress) {
    message("constructing connected subgraphs:")
    fcall = pbapply::pblapply
  } else {
    fcall = lapply
  }

  # determine distinct zones
  czones = unlist(do.call(fcall, fcall_list), 
                  use.names = FALSE, recursive = FALSE)
  czones[distinct(czones)]
}

arg_check_rflex_zones = function(nn, w, cases, ex, 
                                 alpha1, type, pop, 
                                 progress) {
  if (is.null(dim(nn))) stop("nn must be matrix-like with non-NULL dim")
  N = nrow(nn)
  if (nrow(w) != N) stop("nrow(w) must match nrow(nn)")
  if (length(cases) != N) stop("length(cases) must match nrow(nn)")
  if (length(alpha1) != 1 | alpha1 <= 0) {
    stop("alpha1 must be in (0, 1]")
  } 
  if (!is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
  if (!is.null(pop)) {
    if (length(pop) != N) {
      stop("length(pop) != nrow(nn)")
    }
  }
  if (length(progress) != 1 | !is.logical(progress)) {
    stop("progress must be a logical value")
  }
}
