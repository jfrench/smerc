#' Double Connected spatial scan test
#' 
#' \code{dc.test} implements the Double Connected spatial 
#' scan test of Costa et al. (2012). 
#' Starting with a single region as a current zone, new 
#' candidate zones are constructed by combining the current 
#' zone with the connected region that maximizes the 
#' resulting likelihood ratio test statistic, with the 
#' added constraint that the region must have at least
#' two connection (i.e., shares a border with) at least
#' two of the regoins in the current zone.  This  
#' procedure is repeated until adding a connected region
#' does not increase the test statistic (or the population or distance 
#' upper
#' bounds are reached).  The same procedure is repeated for 
#' each region.  The 
#' clusters returned are non-overlapping, ordered from most 
#' significant to least significant.  The first cluster is 
#' the most likely to be a cluster.  If no significant 
#' clusters are found, then the most likely cluster is 
#' returned (along with a warning).
#' 
#' The maximum intercentroid distance can be found by 
#' executing the command: 
#' \code{sp::spDists(as.matrix(coords), longlat = lonlat)}, 
#' based on the specified values of \code{coords} and 
#' \code{lonlat}.
#' 
#' @inheritParams dmst.test
#' @return Returns a list of length two of class scan. The
#'   first element (clusters) is a list containing the
#'   significant, non-ovlappering clusters, and has the the
#'   following components: \item{locids}{The location ids of
#'   regions in a significant cluster.} \item{pop}{The total
#'   population in the cluser window.} \item{cases}{The
#'   observed number of cases in the cluster window.} 
#'   \item{expected}{The expected number of cases in the
#'   cluster window.} \item{smr}{Standarized mortaility
#'   ratio (observed/expected) in the cluster window.} 
#'   \item{rr}{Relative risk in the cluster window.} 
#'   \item{loglikrat}{The loglikelihood ratio for the
#'   cluster window (i.e., the log of the test statistic).} 
#'   \item{pvalue}{The pvalue of the test statistic
#'   associated with the cluster window.} The second element
#'   of the list is the centroid coordinates.  This is
#'   needed for plotting purposes.
#' @author Joshua French
#' @export
#' @seealso \code{\link{scan.stat}},
#'   \code{\link{plot.scan}}, \code{\link{scan.test}},
#'   \code{\link{flex.test}}, \code{\link{uls.test}},
#'   \code{\link{bn.test}}
#' @references Costa, M.A. and Assuncao, R.M. and Kulldorff, M. (2012)
#'   Constrained spanning tree algorithms for
#'   irregularly-shaped spatial clustering, Computational
#'   Statistics & Data Analysis, 56(6), 1771-1783.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' \dontrun{
#' out = dc.test(coords = coords, cases = floor(nydf$cases), 
#'                  pop = nydf$pop, w = nyw, 
#'                  alpha = 0.12, lonlat = TRUE,
#'                  nsim = 5, ubpop = 0.1, ubd = 0.2)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))}
dc.test = function(coords, cases, pop, w,
                   ex = sum(cases)/sum(pop)*pop,
                   nsim = 499, alpha = 0.1, 
                   ubpop = 0.5, ubd = 1,
                   lonlat = FALSE, cl = NULL) {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nsim + 1, ubpop, lonlat, FALSE, 
                      k = 1, w = w)
  
  coords = as.matrix(coords) # ensure proper format
  N = nrow(coords) # number of regions
  ty = sum(cases) # sum of all cases

  # intercentroid distances
  d = sp::spDists(as.matrix(coords), longlat = TRUE)
  # upperbound for distance between centroids in zone
  max_dist = ubd * max(d)
  max_pop = ubpop * sum(pop)
  # find all neighbors from each starting zone within distance upperbound
  all_neighbors = lapply(seq_along(cases), function(i) which(d[i,] <= max_dist))
  
  max_zones = lapply(seq_along(cases), function(i) {
    mst.seq(i, all_neighbors[[i]], cases = cases, 
            pop = pop, w = w, ex = es, ty = ty, 
            max_pop = max_pop, type = "pruned", 
            early = TRUE, nlinks = "two")
  })
  
  # extract statistics from each zone
  tobs = sapply(max_zones, getElement, name = "loglikrat")
  
  if (nsim > 0) {
    # determine which call for simulations
    fcall = pbapply::pblapply
    # setup list for call
    fcall_list = list(X = as.list(seq_len(nsim)), FUN = function(i){
      # simulate new data set
      ysim = stats::rmultinom(1, size = ty, prob = ex)
      # find max statistics for each candidate zone
      tall = mst.all(neighbors = all_neighbors, 
             cases = ysim, pop = pop, w = w, ex = ex,
             ty = ty, max_pop = max_pop, type = "maxonly", 
             early = TRUE, nlinks = "two")
      # return max of statistics for simulation
      return(max(tall))
    }, cl = cl)
    
    # get max statistics for simulated data sets
    tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
    
    # p-values associated with these max statistics for each centroid
    pvalue = unname(sapply(tobs, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  } else {
    pvalue = rep(1, length(tobs))
  }
  
  # determine which potential clusters are significant
  sigc = which(pvalue <= alpha, useNames = FALSE)
  
  # if there are no significant clusters, return most likely cluster
  if (length(sigc) == 0)
  {
    sigc = which.max(tobs)
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # which statistics are significant
  sig_tscan = unlist(tobs, use.names = FALSE)[sigc]
  # order statistics from smallest to largest
  o_sig = order(sig_tscan, decreasing = TRUE)
  # idx of significant clusters in order of significance
  sigc = sigc[o_sig]
  
  # determine the location ids in each significant cluster
  sig_regions = lapply(sigc, function(i) max_zones[[i]]$locids)
  # determine idx of unique non-overlapping clusters
  u = smacpod::noc(sig_regions)
  # return non-overlapping clusters (in order of significance)
  sig_regions = sig_regions[u]
  # unique significant clusters (in order of significance)
  usigc = sigc[u]
  
  # for the unique, non-overlapping clusters in order of significance,
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ration, relative risk,
  sig_tstat = tobs[usigc]
  sig_p = pvalue[usigc]
  sig_coords = coords[usigc,, drop = FALSE]
  sig_r = d[cbind(usigc, sapply(sig_regions, utils::tail, n = 1))]
  sig_yin = sapply(usigc, function(x) max_zones[[x]]$cases)
  sig_ein = sapply(usigc, function(x) max_zones[[x]]$expected)
  sig_popin = sapply(usigc, function(x) max_zones[[x]]$population)
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))

  # reformat output for return
  clusters = vector("list", length(u))
  for (i in seq_along(clusters))
  {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
    clusters[[i]]$r = sig_r[i]
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$loglikrat = sig_tstat[[i]]
    clusters[[i]]$pvalue = sig_p[i]
    clusters[[i]]$w = w[sig_regions[[i]], sig_regions[[i]]]
  }
  
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}


