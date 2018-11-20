#' Determine zones for the Maximum Linking scan test
#' 
#' \code{mlink.zones} determines the zones for the Double 
#' Connection scan test (\code{\link{mlink.test}}).  The 
#' function returns the zones, as well as the associated
#' test statistic, cases in each zone, the expected number
#' of cases in each zone, and the population in each zone.
#' 
#' Every zone considered must have a total population less
#' than \code{ubpop * sum(pop)}.  Additionally, the maximum
#' intercentroid distance for the regions within a zone must
#' be no more than \code{ubd * the maximum intercentroid
#' distance across all regions}.
#' @inheritParams mlink.test
#' @inheritParams mst.all
#' @inheritParams flex.zones
#' @inherit dc.zones return
#' @author Joshua French
#' @references Costa, M.A. and Assuncao, R.M. and Kulldorff, M. (2012)
#'   Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783. 
#'   <https://doi.org/10.1016/j.csda.2011.11.001>
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' # find zone with max statistic starting from each individual region
#' all_zones = mlink.zones(coords, cases = floor(nydf$cases),
#'                         nydf$pop, w = nyw, ubpop = 0.25,
#'                         ubd = .25, longlat = TRUE)
mlink.zones = function(coords, cases, pop, w, 
                       ex = sum(cases)/sum(pop)*pop, 
                       ubpop = 0.5, ubd = 1, longlat = FALSE, 
                       cl = NULL, progress = TRUE) {
  # sanity checking
  arg_check_dmst_zones(coords = coords, cases = cases, 
                       pop = pop, w = w, ex = ex, 
                       ubpop = ubpop, ubd = ubd, 
                       longlat = longlat, type = "all", 
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

  out = mst.all(neighbors = all_neighbors, cases = cases, 
          pop = pop, w = w,  
          ex = ex, ty = ty, max_pop = max_pop, 
          type = "all", nlinks = "max",
          early = FALSE, cl = cl, 
          progress = progress)
  nn = lapply(out, getElement, name = "locids")
  zones = unlist(lapply(nn, function(x) sapply(seq_along(x), function(i) x[seq_len(i)])), recursive = FALSE)
  loglikrat = unlist(lapply(out, getElement, name = "loglikrat"))
  cases = unlist(lapply(out, getElement, name = "cases"))
  expected = unlist(lapply(out, getElement, name = "expected"))
  population = unlist(lapply(out, getElement, name = "population"))
  return(list(zones = zones,
              loglikrat = loglikrat,
              cases = cases,
              expected = expected,
              population = population))
}