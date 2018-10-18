#' Restricted Flexibly-shaped Spatial Scan Test
#' 
#' \code{rflex.test} performs the restricted flexibly shaped spatial scan test of Tango and Takahashi (2012).
#' 
#' The test is performed using the spatial scan test based on the Poisson test statistic and a fixed number of cases.  The first cluster is the most likely to be a cluster.  If no significant clusters are found, then the most likely cluster is returned (along with a warning).
#' 
#' @inheritParams flex.test
#' @inheritParams rflex.zones
#'
#' @return Returns a list of length two of class scan. The first element (clusters) is a list containing the significant, non-ovlappering clusters, and has the the following components: 
#' \item{coords}{The centroid of the significant clusters.}
#' \item{r}{The radius of the window of the clusters.}
#' \item{pop}{The total population in the cluser window.}
#' \item{cases}{The observed number of cases in the cluster window.}
#' \item{expected}{The expected number of cases in the cluster window.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster window.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{loglikrat}{The loglikelihood ratio for the cluster window (i.e., the log of the test statistic).}
#' \item{pvalue}{The pvalue of the test statistic associated with the cluster window.}
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @export
#' @seealso \code{\link{scan.stat}}, \code{\link{plot.scan}}, 
#' \code{\link{scan.test}}, \code{\link{uls.test}}, 
#' \code{\link{dmst.test}}, \code{\link{bn.test}}
#' @references Tango, T., & Takahashi, K. (2005). A flexibly shaped spatial scan statistic for detecting clusters. International journal of health geographics, 4(1), 11.  Kulldorff, M. (1997) A spatial scan statistic. Communications in Statistics -- Theory and Methods 26, 1481-1496.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = rflex.test(coords = coords, cases = floor(nydf$cases),
#'                  w = nyw, k = 10,  
#'                  pop = nydf$pop, nsim = 49, 
#'                  alpha = 0.05, longlat = TRUE)
#'                 
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
rflex.test = function(coords, cases, pop, w, k = 50, 
                     ex = sum(cases)/sum(pop)*pop, 
                     type = "poisson",
                     nsim = 499, 
                     alpha = 0.1, 
                     longlat = FALSE, 
                     alpha1 = 0.2,
                     cl = NULL) {
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nsim + 1, 0.5, longlat, FALSE, k = k, w = w)
  coords = as.matrix(coords)
  N = nrow(coords)

  # compute k nearest neighbors
  nn = knn(coords = coords, longlat = longlat, k = k)  
  
  # determine zones for observed data
  zones = rflex.zones(nn = nn, w = w, 
                      cases = cases, ex = ex, alpha1 = alpha1,
                      cl = cl, progress = FALSE)
  # compute needed information
  ein = unlist(lapply(zones, function(x) sum(ex[x])), use.names = FALSE)
  ty = sum(cases)
  eout = ty - ein
  yin = unlist(lapply(zones, function(x) sum(cases[x])))
  # compute observed scan statistics
  tobs = scan.stat(yin, ein, ty - ein, ty, type = type)
  tscan = max(tobs)
  
  if (nsim > 1) {
    tsim = rflex.sim(nsim = nsim, nn = nn, w = w, ex = ex, 
                     type = type, alpha1 = alpha1, cl = cl)
    pvalue = unname(sapply(tobs, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  } else {
    pvalue = rep(1, length(tobs))
  }
  
  # if (nsim >= 1) {
  #   tsim = pbapply::pbsapply(X = seq_len(nsim),
  #                            FUN = rflex.sim,
  #                            nn = nn, w = w,
  #                            ex = ex, type = type,
  #                            alpha1 = alpha1, cl = cl)
  #   pvalue = unname(sapply(tobs, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  # } else {
  #   pvalue = rep(1, length(tobs))
  # }

  # determine which potential clusters are significant
  sigc = which(pvalue <= alpha, useNames = FALSE)
  
  # if there are no significant clusters, return most likely cluster
  if (length(sigc) == 0) {
    sigc = which.max(tobs)
    warning("No significant clusters.  Returning most likely cluster.")
  }
  
  # only keep significant clusters
  zones = zones[sigc]
  tobs = tobs[sigc]
  
  # order zones from largest to smallest test statistic
  ozones = order(tobs, decreasing = TRUE)
  zones = zones[ozones]
  tobs = tobs[ozones]
  
  # determine significant non-overlapping clusters
  sig = smacpod::noc(zones)
  
  # for the unique, non-overlapping clusters in order of significance,
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ration, relative risk,
  sig_regions = zones[sig]
  sig_tstat = tobs[sig]
  sig_p = pvalue[ozones[sig]]
  centroid = sapply(sig_regions, utils::head, n = 1)
  boundary = sapply(sig_regions, utils::tail, n = 1)
  sig_coords = coords[sapply(sig_regions, function(x) x[1]),, drop = FALSE]
  sig_yin = sapply(sig_regions, function(x) sum(cases[x]))
  sig_ein = sapply(sig_regions, function(x) sum(ex[x]))
  sig_popin = sapply(sig_regions, function(x) sum(pop[x]))
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  sig_w = sapply(centroid, function(x) w[x,,drop = FALSE], simplify = FALSE)

  # reformat output for return
  clusters = vector("list", length(sig_regions))
  for (i in seq_along(clusters)) {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
    clusters[[i]]$r = max(sp::spDists(coords[sig_regions[[i]], ,drop = FALSE]), longlat = longlat)
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$loglikrat = sig_tstat[[i]]
    clusters[[i]]$pvalue = sig_p[i]
    clusters[[i]]$w = sig_w[[i]]
  }
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}

