#' Maxima Likelihood First Scan Test
#' 
#' \code{mlf.test} implements the Maxima Likelihood First scan test of Yao et al. (2011), which is actually a special case of the Dynamic Minimum Spanning Tree of Assuncao et al. (2006).  Find the single region that maximizes the likelihood ratio test statistic.  Starting with this single region as a current zone, new candidate zones are constructed by combining the current zone with the connected region that maximizes the likelihood ratio test static.  This procedure is repeated until the population upper bound is reached.
#' 
#' Only a single cluster is ever returned because the algorithm only constructs a single sequence of starting zones, and overlapping zones are not returned.  Only the zone that maximizes the likelihood ratio test statistic is returned.
#' 
#' @param coords An \eqn{n \times 2} matrix of centroid coordinates for the regions.
#' @param cases The number of cases in each region.
#' @param pop The population size of each region.
#' @param w The binary spatial adjacency matrix.
#' @param ex The expected number of cases for each region.  The default is calculated under the constant risk hypothesis.  
#' @param nsim The number of simulations from which to compute p-value.
#' @param nreport The frequency with which to report simulation progress.  The default is \code{nsim+ 1}, meaning no progress will be displayed.
#' @param ubpop The upperbound of the proportion of the total population to consider for a cluster.
#' @param ubd The upperbound for the proportion of the maximum intercentroid distance to allow for the maximum size of a zone.
#' @param alpha The significance level to determine whether a cluster is signficant.  Default is 0.05.
#' @param lonlat If lonlat is TRUE, then the great circle distance is used to calculate the intercentroid distance.  The default is FALSE, which specifies that Euclidean distance should be used.
#' @param parallel A logical indicating whether the test should be parallelized using the \code{parallel::mclapply function}.  Default is TRUE.  If TRUE, no progress will be reported.
#'
#' @return Returns a list of length two of class scan. The first element (clusters) is a list containing the significant, non-ovlappering clusters, and has the the following components: 
#' \item{locids}{The location ids of regions in a significant cluster.}
#' \item{pop}{The total population in the cluser window.}
#' \item{cases}{The observed number of cases in the cluster window.}
#' \item{expected}{The expected number of cases in the cluster window.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster window.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{loglikrat}{The loglikelihood ratio for the cluster window (i.e., the log of the test statistic).}
#' \item{pvalue}{The pvalue of the test statistic associated with the cluster window.}
#' \item{w}{The adjacency matrix of the cluster.}
#' \item{r}{The maximum radius of the cluster (in terms of intercentroid distance from the starting region).}
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @importFrom parallel mclapply
#' @importFrom stats rmultinom
#' @export
#' @references Yao, Z., Tang, J., & Zhan, F. B. (2011). Detection of arbitrarily-shaped clusters using a neighbor-expanding approach: A case study on murine typhus in South Texas. International journal of health geographics, 10(1), 1.
#' 
#' Assuncao, R.M., Costa, M.A., Tavares, A. and Neto, S.J.F. (2006). Fast detection of arbitrarily shaped disease clusters, Statistics in Medicine, 25, 723-742.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = mlf.test(coords = coords, cases = floor(nydf$cases), 
#'                   pop = nydf$pop, w = nyw, 
#'                   alpha = 0.12, lonlat = TRUE,
#'                   nsim = 10, ubpop = 0.1, ubd = 0.5)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
mlf.test = function (coords, cases, pop, w, 
                     ex = sum(cases)/sum(pop)*pop,
                     nsim = 499, alpha = 0.1, 
                     nreport = nsim + 1, 
                     ubpop = 0.5, ubd = 0.5,
                     lonlat = FALSE, parallel = TRUE) 
{
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nreport, ubpop, lonlat, parallel, 
                      k = 1, w = w)
  
  coords = as.matrix(coords)
  N = nrow(coords)
  ty = sum(cases) # sum of all cases
  
  # calculate test statistics for each individual region
  # yin, ein, eout, ty
  eout = ty - ex
  tobs = scan.stat(cases, ex, eout, ty)
  
  # determine starting region for maxima likelihood first algorithm
  start = which.max(tobs)
  
  # intercentroid distances 
  d = sp::spDists(coords, longlat = lonlat)
  
  # upperbound for population in zone
  max_pop = ubpop *sum(pop)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)
  
  # find neighbors of all regions
  all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
  
  # return sequence of candidate zones (or a subset depending on type)
  max_zone = dmst_max_zone(start, all_neighbors[[start]], cases, pop, w, ex, ty, max_pop, "pruned")
  
  # determine which call for simulations
  fcall = lapply
  if (parallel) fcall = parallel::mclapply
  # setup list for call
  fcall_list = list(X = as.list(seq_len(nsim)), FUN = function(i){
    # simulate new data set
    ysim = stats::rmultinom(1, size = ty, prob = ex)
    # 
    sim_tstat = scan.stat(ysim, ex, eout, ty)
    # determine starting region for maxima likelihood first algorithm
    sim_start = which.max(sim_tstat)

    # update progress
    if ((i%%nreport) == 0) cat(paste(i, ""))
    
    # find max statistic for best candidate zone
    dmst_max_zone(sim_start, all_neighbors[[sim_start]], cases, pop, w, ex, ty, max_pop, type = "maxonly")
  })
  
  # get max statistics for simulated data sets
  tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  
  # p-values associated with these max statistics for each centroid
  pvalue = (sum(tsim >= max_zone$loglikrat) + 1)/(nsim + 1)
  
   # if there are no significant clusters, return most likely cluster
  if(pvalue >= alpha | nsim == 0)
  {
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # for the most likely cluster:
  # find the zone radius, p-value, standarized mortality ratio, relative risk, adjacency matrix
  # and max radius
  startpt = coords[start,] # starting region
  max_zone$r = max(sp::spDistsN1(coords[max_zone$locids,], startpt, longlat = lonlat))
  max_zone$pvalue = pvalue
  max_zone$smr = max_zone$cases/max_zone$expected
  max_zone$rr = (max_zone$cases/max_zone$pop)/((ty - max_zone$cases)/(sum(pop) - max_zone$pop))
  max_zone$w = w[max_zone$locids, max_zone$locids]
  
  # reformat output for return
  outlist = list(clusters = list(max_zone), coords = coords)
  class(outlist) = "scan"
  return(outlist)
}


