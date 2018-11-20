#' Determine zones for the Double Connection scan test
#' 
#' \code{dc.zones} determines the zones for the Double 
#' Connection scan test (\code{\link{dc.test}}).
#' 
#' Every zone considered must have a total population less
#' than \code{ubpop * sum(pop)}.  Additionally, the maximum
#' intercentroid distance for the regions within a zone must
#' be no more than \code{ubd * the maximum intercentroid
#' distance across all regions}.
#' 
#' \code{type} is a character vector indicating what should 
#' be returned by the function.  If \code{type = "maxonly"},
#' then the maximum test statistic from each starting region
#' is returned .  If \code{type = "pruned"}, the function
#' returns a list that includes the location ids, test
#' statistic, total cases, expected cases, and total
#' population for the zone with the maximum test statistic
#' for each starting region.  If \code{type = "all"}, the
#' function returns a list of lists that includes the
#' location ids, test statistic, total cases, expected
#' cases, and total population for the sequence of candidate
#' zones associated with each starting region.
#' @inheritParams dmst.test
#' @inheritParams mst.all
#' @inheritParams flex.zones
#' @return Returns a list of relevant information.  See 
#'   Details.
#' @author Joshua French
#' @references Costa, M.A. and Assuncao, R.M. and Kulldorff, M. (2012)
#'   Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783. 
#'   <https://doi.org/10.1016/j.csda.2011.11.001>
#'   
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' # find zone with max statistic starting from each individual region
#' max_zones = dc.zones(coords, cases = floor(nydf$cases),
#'                      nydf$pop, w = nyw, ubpop = 0.25,
#'                      ubd = .25, longlat = TRUE)
dc.zones = function(coords, cases, pop, w, 
                    ex = sum(cases)/sum(pop)*pop, 
                    ubpop = 0.5, ubd = 1, longlat = FALSE, 
                    type = "maxonly", cl = NULL,
                    progress = TRUE) {
  # sanity checking
  arg_check_dmst_zones(coords = coords, cases = cases, 
                       pop = pop, w = w, ex = ex, 
                       ubpop = ubpop, ubd = ubd, 
                       longlat = longlat, type = type, 
                       progress = progress)  
  # setup various arguments and such
  ty = sum(cases)   # total number of cases
  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  # upperbound for population in zone
  max_pop = ubpop * sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)
  # find all neighbors from each starting zone within distance upperbound
  all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))

  mst.all(neighbors = all_neighbors, cases = cases, 
          pop = pop, w = w,  
          ex = ex, ty = ty, max_pop = max_pop, 
          type = type, nlinks = "two",
          early = TRUE, cl = cl, 
          progress = progress)
}