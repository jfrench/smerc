#' Upper Level Set Spatial Scan Test
#'
#' \code{uls.test} performs the Upper Level Set (ULS)
#' spatial scan test of Patil and Taillie (2004).  The test
#' is performed using the spatial scan test based on a fixed
#' number of cases.  The windows are based on the Upper
#' Level Sets proposed by Patil and Taillie (2004).  The
#' clusters returned are non-overlapping, ordered from most
#' significant to least significant.  The first cluster is
#' the most likely to be a cluster.  If no significant
#' clusters are found, then the most likely cluster is
#' returned (along with a warning).
#'
#' The ULS method has a special (and time consuming)
#' construction when the observed rates aren't unique.  This
#' is unlikely to arise for real data, except with observed
#' rates of 0, which are of little interest.  The method can
#' take substantially if this is considered.
#'
#' @param w A binary spatial adjacency matrix for the
#'   regions.
#' @param check.unique A logical value indicating whether a
#'   check for unique values should be determined.  The
#'   default is \code{FALSE}.  This is unlikely to make a
#'   practical different for most real data sets.
#' @inheritParams scan.test
#'
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
#' @seealso \code{\link{scan.stat}},
#'   \code{\link{plot.scan}}, \code{\link{scan.test}},
#'   \code{\link{flex.test}}, \code{\link{dmst.test}},
#'   \code{\link{bn.test}}
#' @export
#' @references Patil, G.P. & Taillie, C. Upper level set
#'   scan statistic for detecting arbitrarily shaped
#'   hotspots. Environmental and Ecological Statistics
#'   (2004) 11(2):183-197.
#'   <doi:10.1023/B:EEST.0000027208.48919.7e>
#' @examples
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = uls.test(coords = coords, cases = floor(nydf$cases),
#'                pop = nydf$pop, w = nyw,
#'                alpha = 0.05, longlat = TRUE,
#'                nsim = 9, ubpop = 0.5)
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
uls.test = function(coords, cases, pop, w,
                    ex = sum(cases)/sum(pop)*pop,  
                    nsim = 499, alpha = 0.1, 
                    ubpop = 0.5, longlat = FALSE,
                    cl = NULL, type = "poisson",
                    check.unique = FALSE) {
  # sanity checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nsim + 1, ubpop, longlat, TRUE, 
                      k = 1, w = w, type = type)
  
  coords = as.matrix(coords)
  zones = uls.zones(cases, pop, w, ubpop)
  
  # compute needed information
  ty = sum(cases)
  yin = zones.sum(zones, cases)
  
  # compute test statistics for observed data
  if (type == "poisson") {
    ein = zones.sum(zones, ex)
    tobs = stat.poisson(yin, ty - yin, ein, ty - ein)
  } else if (type == "binomial") {
    tpop = sum(pop)
    popin = zones.sum(zones, pop)
    tobs = stat.binom(yin, ty - yin, ty, popin, tpop - popin, tpop)
  }
  
  # compute test statistics for simulated data
  if (nsim > 0) {
    message("computing statistics for simulated data:")
    tsim = uls.sim(nsim = nsim, ty = ty, ex = ex, w = w,
                   pop = pop, ubpop = ubpop, cl = cl)
    pvalue = mc.pvalue(tobs, tsim)
  } else {
    pvalue = rep(1, length(tobs))
  }
  
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
  pvalue = pvalue[sigc]
  
  prep.scan(tobs = tobs, zones = zones, pvalue = pvalue, 
            coords = coords, cases = cases, pop = pop,
            ex = ex, longlat = longlat, w = w,
            d = NULL)
  # coords = as.matrix(coords)
  # N = nrow(coords)
  # 
  # ty = sum(y) # sum of all cases
  # 
  # fcall = pbapply::pblapply
  # fcall_list = list(X = as.list(1:nsim), FUN = function(i){
  #   # simulate new data set
  #   ysim = stats::rmultinom(1, size = ty, prob = e)
  #   uz = uls.zones(ysim, pop, w, ubpop = ubpop)
  #   # cumulate the number of cases inside the successive windows
  #   yin = unlist(lapply(uz, function(x) sum(ysim[x])), use.names = FALSE)
  #   ein = unlist(lapply(uz, function(x) sum(e[x])), use.names = FALSE)
  #   # calculate all test statistics
  #   tall = scan.stat(yin, ein, ty - ein, ty)
  #   # return max of statistics for simulation
  #   return(max(tall))
  # }, cl = cl)
  # 
  # # use mclapply or lapply to find max statistics for each simulation
  # tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
  # 
  # # determine yin and yout for all windows for observed data
  # uz = uls.zones(y, pop, w, ubpop)
  # yin = unlist(lapply(uz, function(x) sum(y[x])), use.names = FALSE)
  # yout = ty - yin
  # ein = unlist(lapply(uz, function(x) sum(e[x])), use.names = FALSE)
  # eout = ty - ein
  # popin = unlist(lapply(uz, function(x) sum(pop[x])), use.names = FALSE)
  # 
  # ### calculate scan statistics for observed data
  # # of distance from observation centroid
  # tobs = scan.stat(yin, ein, ty - ein, ty)
  # # max scan statistic over all windows
  # tscan = max(tobs)
  # 
  # # p-values associated with these max statistics for each centroid
  # pvalue = unname(sapply(tobs, function(x) (sum(tsim >= x) + 1)/(nsim + 1)))
  # 
  # # determine which potential clusters are significant
  # sigc = which(pvalue <= alpha, useNames = FALSE)
  # 
  # # if there are no significant clusters, return most likely cluster
  # if (length(sigc) == 0) {
  #   sigc = which.max(tobs)
  #   warning("No significant clusters.  Returning most likely cluster.")
  # }
  # 
  # # which statistics are significant
  # sig_tscan = tobs[sigc]
  # # order statistics from smallest to largest
  # o_sig = order(sig_tscan, decreasing = TRUE)
  # # idx of significant clusters in order of significance
  # sigc = sigc[o_sig]
  # 
  # # determine the location ids in each significant cluster
  # # sig_regions = mapply(function(a, b) mynn[[a]][1:b], a = sigc, b = tmax_pos[sigc], SIMPLIFY = FALSE) 
  # sig_regions = uz[sigc]
  # # determine idx of unique non-overlapping clusters
  # u = smacpod::noc(sig_regions)
  # # return non-overlapping clusters (in order of significance)
  # sig_regions = sig_regions[u]
  # # unique significant clusters (in order of significance)
  # usigc = sigc[u]
  # 
  # # for the unique, non-overlapping clusters in order of significance,
  # # find the associated test statistic, p-value, centroid,
  # # window radius, cases in window, expected cases in window, 
  # # population in window, standarized mortality ration, relative risk,
  # sig_tstat = tobs[usigc]
  # sig_p = pvalue[usigc]
  # sig_yin = yin[usigc]
  # sig_ein = ein[usigc]
  # sig_popin = popin[usigc]
  # sig_smr = sig_yin/sig_ein
  # sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  # 
  # # reformat output for return
  # clusters = vector("list", length(u))
  # for(i in seq_along(clusters))
  # {
  #   clusters[[i]]$locids = sig_regions[[i]]
  #   clusters[[i]]$pop = sig_popin[i]
  #   clusters[[i]]$cases = sig_yin[i]
  #   clusters[[i]]$expected = sig_ein[i]
  #   clusters[[i]]$smr = sig_smr[i]
  #   clusters[[i]]$rr = sig_rr[i]
  #   clusters[[i]]$loglikrat = sig_tstat[i]
  #   clusters[[i]]$pvalue = sig_p[i]
  #   clusters[[i]]$w = w[sig_regions[[i]], sig_regions[[i]]]
  # }
  # outlist = list(clusters = clusters, coords = coords)
  # class(outlist) = "scan"
  # return(outlist)
}


