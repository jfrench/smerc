#' Elliptical Spatial Scan Test
#'
#' \code{elliptic.test} performs the elliptical scan test of
#' Kulldorf et al. (2006).
#'
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  Candidate zones are elliptical and extend from the
#' observed data locations.  The clusters returned are
#' non-overlapping, ordered from most significant to least
#' significant.  The first cluster is the most likely to
#' be a cluster.  If no significant clusters are found, then
#' the most likely cluster is returned (along with a
#' warning).
#' @inheritParams scan.test
#' @param shape The ratios of the major and minor axes of
#'   the desired ellipses.
#' @param nangle The number of angles (between 0 and 180) to
#'   consider for each shape.
#' @param a The penalty for the spatial scan statistic.  The
#'   default is 0.5.
#' @inheritParams pbapply::pblapply
#'
#' @return Returns a list of length two of class scan. The
#'   first element (clusters) is a list containing the
#'   significant, non-overlappering clusters, and has the
#'   the following components: \item{locids}{The location
#'   ids of regions in a significant cluster.}
#'   \item{coords}{The centroid of the significant
#'   clusters.} \item{r}{The radius of the cluster (the
#'   largest intercentroid distance for regions in the
#'   cluster).} \item{pop}{The total population of the
#'   regions in the cluster.} \item{cases}{The observed
#'   number of cases in the cluster.} \item{expected}{The
#'   expected number of cases in the cluster.}
#'   \item{smr}{Standarized mortaility ratio
#'   (observed/expected) in the cluster.} \item{rr}{Relative
#'   risk in the cluster.} \item{loglikrat}{The
#'   loglikelihood ratio for the cluster (i.e., the log of
#'   the test statistic).} \item{pvalue}{The pvalue of the
#'   test statistic associated with the cluster.} The second
#'   element of the list is the centroid coordinates.  This
#'   is needed for plotting purposes.
#' @seealso \code{\link{scan.stat}}, \code{\link{plot.scan}}
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial scan
#' statistic. Communications in Statistics -- Theory and
#' Methods, 26:1481-1496. 
#' 
#' Kulldorff, M., Huang, L., Pickle,
#' L. and Duczmal, L. (2006) An elliptic spatial scan
#' statistic. Statististics in Medicine, 25:3929-3943.
#' doi:10.1002/sim.2490
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = elliptic.test(coords = coords,
#'                    cases = floor(nydf$cases),
#'                    pop = nydf$pop, ubpop = 0.1, 
#'                    nsim = 49,
#'                    alpha = 0.12,
#'                    shape = 1.5, nangle = 4)
elliptic.test = function(coords, cases, pop, 
                     ex = sum(cases)/sum(pop)*pop, 
                     nsim = 499, alpha = 0.1,  
                     ubpop = 0.5,  
                     shape = c(1, 1.5, 2, 3, 4, 5),
                     nangle = c(1, 4, 6, 9, 12, 15),
                     a = 0.5,
                     cl = NULL,
                     type = "poisson", 
                     min.cases = 2) {
  # argument checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha,
                      nsim + 1, ubpop, TRUE, TRUE,
                      k = 1, w = diag(nrow(coords)))
  if (length(shape) != length(nangle)) {
    stop("The length of shape and nangle must match.")
  }
  
  # convert to proper format
  coords = as.matrix(coords)
  N = nrow(coords)
  # short names

  enn = elliptic.nn(coords, pop = pop, ubpop = ubpop,
                    shape = shape, nangle = nangle)
  # determine the expected cases in/out each successive 
  # window, total number of cases, total population
  ein = unlist(lapply(enn$nn, function(x) cumsum(ex[x])))
  ty = sum(cases) # sum of all cases
  eout = ty - ein # counts expected outside the window

  # determine yin and yout for all windows for observed data
  yin = unlist(lapply(enn$nn, function(x) cumsum(cases[x])))

  ### calculate scan statistics for observed data
  # of distance from observation centroid
  tobs = stat.poisson(yin, ty - yin, ein, eout, a = a, 
                      shape = enn$shape_all)

  # determine distinct zones
  pri = randtoolbox::get.primes(N)
  wdup = duplicated(unlist(lapply(enn$nn, function(x) cumsum(log(pri[x])))))
  
  # determine positions in nn of all zones
  nnn = unlist(lapply(enn$nn, length), use.names = FALSE)
  allpos = cbind(rep(seq_along(enn$nn), times = nnn),
                 unlist(sapply(nnn, seq_len)))
 
  # remove zones with a test statistic of 0
  # or fewer than minimum number of cases or
  # indistinct
  w0 = which(tobs == 0 | yin < min.cases | wdup)
  tobs = tobs[-w0]
  shape_all = enn$shape_all[-w0]
  angle_all = enn$angle_all[-w0]
  allpos = allpos[-w0,]

  # setup list for call
  if (nsim > 0) {
  #   fcall = pbapply::pblapply
  #   fcall_list = list(X = seq_len(nsim), FUN = function(i){
  #     # simulate new data set
  #     ysim = stats::rmultinom(1, size = ty, prob = e)
  #     # cumulate the number of cases inside the successive windows
  #     yin = unlist(lapply(mynn, function(x) cumsum(ysim[x])))
  #     # calculate all test statistics
  #     tall = scan.stat(yin, ein, eout, ty, type, a = a, shape = shape)
  #     # return max of statistics for simulation
  #     return(max(tall))
  #   }, cl = cl)
    
    # # use mclapply or lapply to find max statistics for each simulation
    # tsim = unlist(do.call(fcall, fcall_list), use.names = FALSE)
    tsim = elliptic.sim(nsim = nsim, nn = enn$nn, ty = ty, 
                        ex = ex, a = a, 
                        shape_all = enn$shape_all,
                        ein = ein, eout = eout, cl = cl)
    
    # p-values associated with these max statistics for each centroid
    pvalue = mc.pvalue(tobs, tsim)
  } else {
    pvalue = rep(1, length(tobs))
  }
  
  zones = apply(allpos, 1, function(x) enn$nn[[x[1]]][seq_len(x[2])])
  
  prep.scan(tobs = tobs, zones = zones, pvalue = pvalue, 
            coords = coords, cases = cases, pop = pop,
            ex = ex, longlat = FALSE, w = NULL,
            d = NULL, a = a, shape_all = shape_all,
            angle_all = angle_all)
  
#   # determine which potential clusters are significant
#   sigc = which(pvalue <= alpha, useNames = FALSE)
# 
#   # if there are no significant clusters, return most likely cluster
#   if (length(sigc) == 0) {
#     sigc = which.max(tobs)
#     warning("No significant clusters.  Returning most likely cluster.")
#   }
# 
#   allpos = allpos[sigc, ]
#   tobs = tobs[sigc]
#   shape_all = shape_all[sigc]
#   angle_all = angle_all[sigc]
# 
#   # construct all zones
#   zones = apply(allpos, 1, function(x) enn$nn[[x[1]]][seq_len(x[2])])
# 
#   # order zones by most significant
#   ozones = order(tobs, decreasing = TRUE)
#   zones = zones[ozones]
# 
#   # determine significant non-overlapping clusters
#   sig = smacpod::noc(zones)
#   
#   # for the unique, non-overlapping clusters in order of significance,
#   # find the associated test statistic, p-value, centroid,
#   # window radius, cases in window, expected cases in window, 
#   # population in window, standarized mortality ratio, relative risk,
#   sig_regions = zones[sig]
#   sig_tstat = tobs[ozones[sig]]
#   sig_p = pvalue[ozones[sig]]
#   sig_shape = shape_all[ozones[sig]]
#   sig_angle = angle_all[ozones[sig]]
#   centroid = sapply(sig_regions, utils::head, n = 1)
#   boundary = sapply(sig_regions, utils::tail, n = 1)
#   sig_coords = coords[sapply(sig_regions, function(x) x[1]),, drop = FALSE]
#   sig_yin = sapply(sig_regions, function(x) sum(y[x]))
#   sig_ein = sapply(sig_regions, function(x) sum(e[x]))
#   sig_popin = sapply(sig_regions, function(x) sum(pop[x]))
#   sig_smr = sig_yin/sig_ein
#   sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
#   sig_w = lapply(sig_regions, function(x) {
#     matrix(c(0, rep(1, length(x) - 1)), nrow = 1)  
#   })
#   sig_minor = unname(sapply(seq_along(sig_regions), function(i) {
#     first = sig_regions[[i]][1]
#     last = utils::tail(sig_regions[[i]], 1)
#     dist.ellipse(coords[c(first, last),,drop = FALSE],
#                  shape = sig_shape[i],
#                  angle = sig_angle[i])[1,2]
#   }))
#   sig_major = sig_minor * sig_shape
#   
#   # reformat output for return
#   clusters = vector("list", length(sig_regions))
#   for (i in seq_along(clusters)) {
#     clusters[[i]]$locids = sig_regions[[i]]
#     clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
#     clusters[[i]]$semiminor.axis = sig_minor[i]
#     clusters[[i]]$semimajor.axis = sig_major[i]
#     clusters[[i]]$angle = sig_angle[i]
#     clusters[[i]]$shape = sig_shape[i]
#     clusters[[i]]$pop = sig_popin[i]
#     clusters[[i]]$cases = sig_yin[i]
#     clusters[[i]]$expected = sig_ein[i]
# #    clusters[[i]]$cases.per.100000 = sig_yin[i]/sig_popin[i] * 1e5
#     clusters[[i]]$smr = sig_smr[i]
#     clusters[[i]]$rr = sig_rr[i]
#     clusters[[i]]$test.statistic = sig_tstat[[i]]
#     clusters[[i]]$loglikrat = scan.stat(sig_yin[i], sig_ein[i], ty - sig_ein[i], ty, shape = sig_shape[i])
#     clusters[[i]]$pvalue = sig_p[i]
#     clusters[[i]]$w = sig_w[[i]]
#   }
#   outlist = list(clusters = clusters, coords = coords)
#   class(outlist) = "scan"
#   return(outlist)
}