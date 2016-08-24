#' Static minimum spanning tree scan test
#' 
#' \code{smst.test} implements the Static Minimum Spanning Tree scan test of Assuncao et al. (2006).  Starting with a single region as a current zone, new candidate zones are constructed by combining the current zone with the connected zone that maximizes the likelihood function.  This is procedure repeated until the population upper bound is reached.  The same procedure is repeated for each zone.  The maxima likelihood first scan test proposed by Yao et al. (2011) is an independent variant of this, but only searches from the single region that maximizes the likelihood ratio scan statistic. 
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
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @importFrom parallel mclapply
#' @importFrom smacpod noc
#' @importFrom stats rmultinom
#' @export
#' @references Assuncao, R.M., Costa, M.A., Tavares, A. and Neto, S.J.F. (2006). Fast detection of arbitrarily shaped disease clusters, Statistics in Medicine, 25, 723-742.
#' 
#' Yao, Z., Tang, J., & Zhan, F. B. (2011). Detection of arbitrarily-shaped clusters using a neighbor-expanding approach: A case study on murine typhus in South Texas. International journal of health geographics, 10(1), 1.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = smst.test(coords = coords, cases = floor(nydf$cases), 
#'                   pop = nydf$pop, w = nyw, 
#'                   alpha = 0.12, lonlat = TRUE,
#'                   nsim = 10, ubpop = 0.1, ubd = 0.5)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
smst.test = function (coords, cases, pop, w, 
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
  
  # smst.zones output (czo)
  czo = smst.zones(coords = coords, cases = cases,
                 pop = pop, w = w, ex = ex, ubpop = ubpop, 
                 ubd = ubd, lonlat = lonlat, parallel = parallel,
                 nsim = nsim)

  uz = czo$uz
  tsim = czo$tsim
  
  ### calculate scan statistics for observed data
  # of distance from observation centroid
  tobs = czo$loglikrat
 
  # p-values associated with these max statistics for each centroid
  pvalue = (sum(tsim >= tobs) + 1)/(nsim + 1)
  
   # if there are no significant clusters, return most likely cluster
  if(pvalue >= alpha | nsim == 0)
  {
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # determine the location ids in each significant cluster
  sig_regions = uz

  # for the cluster:
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ration, relative risk,
  sig_tstat = tobs
  sig_p = pvalue
  sig_yin = czo$cases
  sig_ein = czo$expected
  sig_popin = czo$population
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  
  # reformat output for return
  clusters = vector("list", 1)
  for(i in seq_along(clusters))
  {
    clusters[[i]]$locids = sig_regions
    clusters[[i]]$pop = sig_popin
    clusters[[i]]$cases = sig_yin
    clusters[[i]]$expected = sig_ein
    clusters[[i]]$smr = sig_smr
    clusters[[i]]$rr = sig_rr
    clusters[[i]]$loglikrat = sig_tstat
    clusters[[i]]$pvalue = sig_p[i]
    clusters[[i]]$w = w[sig_regions, sig_regions]
  }
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}


