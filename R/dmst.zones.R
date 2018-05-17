#' Determine zones for the dynamic minimum spanning tree
#' scan test of Assuncao et al. (2006)
#' 
#' \code{dmst.zones} determines the zones that produce the
#' largest test statistic using a greedy algorithm. 
#' Specifically, starting individually with each region as a
#' starting zone, new (connected) regions are added to the
#' current zone in the order that results in the largest
#' likelihood ratio test statistic.  This is used to
#' implement the dynamic minimum spanning tree (dmst) scan
#' test of Assuncao et al. (2006).
#' 
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  The first cluster is the most likely to be a
#' cluster.  If no significant clusters are found, then the
#' most likely cluster is returned (along with a warning).
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
#' @return Returns a list of relevant information.  See 
#'   Details.
#' @author Joshua French
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and
#'   Neto, S.J.F. (2006). Fast detection of arbitrarily
#'   shaped disease clusters, Statistics in Medicine, 25,
#'   723-742.
#'   
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' # find zone with max statistic starting from each individual region
#' max_zones = dmst.zones(coords, cases = floor(nydf$cases),
#'                        nydf$pop, w = nyw, ubpop = 0.25,
#'                        ubd = .25, longlat = TRUE)
#' head(max_zones)

# data(nydf)
# data(nyw)
# coords = nydf[,c("longitude", "latitude")]
# cases = floor(nydf$cases)
# pop = nydf$population
# w = nyw
# e = sum(cases)/sum(pop)*pop
# ubpop = 0.5
# ubd = 0.5
# longlat = TRUE
dmst.zones = function(coords, cases, pop, w, 
                      ex = sum(cases)/sum(pop)*pop, 
                      ubpop = 0.5, ubd = 1, longlat = FALSE, 
                      type = "maxonly", 
                      cl = NULL)
{
  # sanity checking
  arg_check_dmst_zones(coords, cases, pop, w, ex, ubpop, ubd, longlat, type)
  
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
  # should only max be returned, or a pruned version

  # set up finding zones with max stat from each starting region
  fcall = pbapply::pblapply
  fcall_list = list(X = as.list(seq_along(all_neighbors)), FUN = function(i){
    # find zone with max stat, starting from each region
    mst.seq(i, all_neighbors[[i]], cases, pop, w, ex, ty, max_pop, type = type)
  }, cl = cl)

  # obtain list of zones with maximum statistic from each starting region (or the max statistic)
  do.call(fcall, fcall_list)
}

arg_check_dmst_zones = function(coords, cases, pop, w, ex, ubpop, ubd, longlat, type)
{
  if(ncol(coords) != 2) stop("coords should have 2 columns")
  if(nrow(coords) != length(cases)) stop("nrow(coords) != length(cases)")
  if(length(cases) != length(pop)) stop('length(cases) != length(pop)')
  if(length(cases) != nrow(w)) stop('length(cases) != nrow(w)')
  if(length(cases) != length(ex)) stop('length(cases) != length(ex)')
  if(length(ubpop) != 1 | !is.numeric(ubpop)) stop("ubpop should be a single number")
  if(ubpop <= 0 | ubpop > 1) stop("ubpop not in (0, 1]")
  if(length(ubd) != 1 | !is.numeric(ubd)) stop("ubd should be a single number")
  if(ubd <= 0 | ubd > 1) stop("ubd not in (0, 1]")
  if(length(longlat) != 1 || !is.logical(longlat)) stop("longlat must be a single logical value")
  if(!is.element(type, c("maxonly", "pruned", "all"))) stop("Invalid type.")
}