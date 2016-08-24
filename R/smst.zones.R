#' Determine zones using the static minimum spanning tree scan test of Assuncao et al. (2006) 
#' 
#' \code{smst.zones} determines the most likely cluster zone obtained by implementing the static minimum spanning tree (SMST) scan test of Assuncao et al. (2006) 
#' 
#' Each step of the SMST scan test seeks to maximize the likelihood ratio test statistic used in the original spatial scan test (Kulldorff 1995).  The first zone considered is the region that maximizes this statistic, providing that no more than \code{ubpop} proportion of the total population is in the zone.  The second zone is the first zone and the connected region that maximizes the scan statistic, subject to the population and distance constraints.
#' 
#' Every zone considered must have a total population less than \code{ubpop * the total population across all regions} in the study area.  Additionally, the maximum intercentroid distance for the regions within a zone must be no more than \code{ubd * the maximum intercentroid distance across all regions}. 
#' 
#' See also the maxima likelihood first scan method of Yao et al. (2011).
#' 
#' @param coords	\eqn{An n \times 2} matrix of centroid coordinates for the regions.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param ex The expected number of cases for each region. The default is calculated under the constant risk hypothesis.
#' @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
#' @param ubd The upperbound of the proportion of the maximum intercentroid distance across all regions.  Default is 1.
#' @param lonlat A logical value indicating whether Euclidean (\code{lonlat = FALSE}) or great circle distance (\code{lonlat = TRUE}) should be used to calculate intercentroid distance.  The default is \code{FALSE}.
#' @param parallel A logical value indicating whether finding candidate zones should use multiple cores.  The default is \code{FALSE}.  Multiple cores are used if \code{paralle = TRUE}.
#' @param maxonly A logical value indicating whether to return only the maximum test statistic across all candidate zones.  Default is \code{FALSE}.

#' @return Returns a list of zones to consider for clustering that includes the location id of each zone and the associated test statistic, counts, expected counts, and population in the zone. If \code{maxonly = TRUE}, then only the maximum test statistic across these zones is returned.
#' @author Joshua French
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and Neto, S.J.F. (2006). Fast detection of arbitrarily shaped disease clusters, Statistics in Medicine, 25, 723-742.
#' 
#' Yao, Z., Tang, J., & Zhan, F. B. (2011). Detection of arbitrarily-shaped clusters using a neighbor-expanding approach: A case study on murine typhus in South Texas. International Journal of Health Geographics, 10(1), 1.
#' @export
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = as.matrix(nydf[,c("x", "y")])
#' smst.zones(coords, cases = floor(nydf$cases), nydf$pop, w = nyw, ubpop = 0.1, ubd = .5)

data(nydf)
data(nyw)
coords = nydf[,c("longitude", "latitude")]
cases = floor(nydf$cases)
pop = nydf$population
w = nyw
ex = sum(cases)/sum(pop)*pop
ubpop = 0.5
ubd = 0.5
lonlat = TRUE
parallel = FALSE

smst.zones = function(coords, cases, pop, w, ex = sum(cases)/sum(pop)*pop, ubpop = 0.5, ubd = 1, lonlat = FALSE, parallel = FALSE)
{
  arg_check_smst_zones(coords, cases, pop, w, ex, ubpop, ubd, lonlat, parallel)
  
  ty = sum(cases)   # total number of cases
  
  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = lonlat)
  
  # upperbound for population in zone
  max_pop = ubpop *sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)
  
  all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
  
  # regions_used = which(pop > max_pop) # exclude certain regions from the get go
  
  # should only max be returned, or a pruned version
  type = ifelse(maxonly, "maxonly", "pruned")
  
  fcall = lapply
  if (parallel) fcall = parallel::mclapply
  fcall_list = list(X = as.list(seq_len(nsim)), FUN = function(i){
    # ysim = stats::rmultinom(1, size = ty, prob = ex)
    # return max statistic for mlf zone
    
    
    smst_max_zone = function(neighbors, cases, pop, w, e, ty, max_pop, type = "maxonly")
    
    smst_max_zone(all_neighbors[[i]], cases, pop, w, e, ty, max_pop, maxonly = TRUE)
    
  })
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  
  # calculate regionwise test statistics for observed data
  stat = scan.stat(cases, ex, ty - ex, ty) * (pop <= max_pop)
  # find starting zone (the region with the largest stat)
  which_max_stat = which.max(stat)
  # vector of unique zones, potential neighbors for unique zones based on distance
  uz =  max_neighbors = vector("list", length(all_neighbors[[which_max_stat]]))
  # first zone
  uz[[1]] = which_max_stat
  # update various quantities
  loglikrat = yin = ein = popin = numeric(length(uz))
  loglikrat[1] = stat[uz[[1]]]
  yin[1] = cases[uz[[1]]]
  ein[1] = ex[uz[[1]]]
  popin[1] = pop[uz[[1]]]
  max_neighbors[[1]] = setdiff(all_neighbors[[uz[[1]]]], c(regions_used, uz[[1]]))
  
  # return most likely zone from smst search
  smst_zone = smst_zones_internal(uz, max_neighbors, cases, pop, w, ex, ty, max_pop, loglikrat, yin, ein, popin, maxonly = FALSE)
  mlf_zone$tsim = tsim
  return(mlf_zone)
}


# smst.zones = function(coords, cases, pop, w, ex = sum(cases)/sum(pop)*pop, ubpop = 0.5, ubd = 1, lonlat = FALSE, parallel = TRUE, nsim = 0)
# {
#   arg_check_mlf_zones(coords, cases, pop, w, ex, ubpop, ubd, lonlat, parallel, nsim)
#   
#   ty = sum(cases)   # total number of cases
# 
#   # intercentroid distances
#   d = sp::spDists(as.matrix(coords), longlat = lonlat)
# 
#   # upperbound for population in zone
#   max_pop = ubpop *sum(pop)
#   # upperbound for distance between centroids in zone
#   max_dist = ubd * max(d)
# 
#   all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
#   
#   regions_used = which(pop > max_pop) # exclude certain regions from the get go
#   
#   fcall = lapply
#   if (parallel) fcall = parallel::mclapply
#   fcall_list = list(X = as.list(seq_len(nsim)), FUN = function(i){
#     ysim = stats::rmultinom(1, size = ty, prob = ex)
#     # scan statistics for each individual region
#     stat = scan.stat(ysim, ex, ty - ex, ty) * (pop <= max_pop)
#     # find starting zone (the region with the largest stat)
#     which_max_stat = which.max(stat)
#     # vector of unique zones, potential neighbors for unique zones based on distance
#     uz =  max_neighbors = vector("list", length(all_neighbors[[which_max_stat]]))
#     # first zone
#     uz[[1]] = which_max_stat
#     # update various quantities
#     loglikrat = yin = ein = popin = numeric(length(uz))
#     loglikrat[1] = stat[uz[[1]]]
#     yin[1] = cases[uz[[1]]]
#     ein[1] = ex[uz[[1]]]
#     popin[1] = pop[uz[[1]]]
#     max_neighbors[[1]] = setdiff(all_neighbors[[uz[[1]]]], c(regions_used, uz[[1]]))
#     
#     # return max statistic for mlf zone
#     mlf_zones_internal(uz, max_neighbors, cases, pop, w, ex, ty, max_pop, loglikrat, yin, ein, popin)
#   })
#   tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
#   
#   # calculate regionwise test statistics for observed data
#   stat = scan.stat(cases, ex, ty - ex, ty) * (pop <= max_pop)
#   # find starting zone (the region with the largest stat)
#   which_max_stat = which.max(stat)
#   # vector of unique zones, potential neighbors for unique zones based on distance
#   uz =  max_neighbors = vector("list", length(all_neighbors[[which_max_stat]]))
#   # first zone
#   uz[[1]] = which_max_stat
#   # update various quantities
#   loglikrat = yin = ein = popin = numeric(length(uz))
#   loglikrat[1] = stat[uz[[1]]]
#   yin[1] = cases[uz[[1]]]
#   ein[1] = ex[uz[[1]]]
#   popin[1] = pop[uz[[1]]]
#   max_neighbors[[1]] = setdiff(all_neighbors[[uz[[1]]]], c(regions_used, uz[[1]]))
#   
#   # return most likely zone from mlf search
#   mlf_zone = mlf_zones_internal(uz, max_neighbors, cases, pop, w, ex, ty, max_pop, loglikrat, yin, ein, popin, maxonly = FALSE)
#   mlf_zone$tsim = tsim
#   return(mlf_zone)
# }

arg_check_smst_zones = function(coords, cases, pop, w, ex, ubpop, ubd, lonlat)
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
  if(length(lonlat) != 1 || !is.logical(lonlat)) stop("lonlat must be a single logical value")
  if(length(parallel) != 1 || !is.logical(parallel)) stop("parallel must be a single logical value")
}