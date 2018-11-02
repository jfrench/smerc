#' Spatial Scan Test
#'
#' \code{scan.test} performs the spatial scan test of
#' Kulldorf (1997).
#'
#' The test is performed using the spatial scan test based
#' on the Poisson test statistic and a fixed number of
#' cases.  Candidate zones are circular and extend from the
#' observed data locations.  The clusters returned are
#' non-overlapping, ordered from most significant to least
#' significant.  The first cluster is the most  cflikely to
#' be a cluster.  If no significant clusters are found, then
#' the most likely cluster is returned (along with a
#' warning).
#'
#' @param coords An \eqn{n \times 2} matrix of centroid
#'   coordinates for the regions.
#' @param cases The number of cases observed in each region.
#' @param pop The population size associated with each
#'   region.
#' @param ex The expected number of cases for each region.
#'   The default is calculated under the constant risk
#'   hypothesis.
#' @param nsim The number of simulations from which to
#'   compute the p-value.
#' @param ubpop The upperbound of the proportion of the
#'   total population to consider for a cluster.
#' @param alpha The significance level to determine whether
#'   a cluster is signficant.  Default is \code{0.10}.
#' @param longlat The default is \code{FALSE}, which
#'   specifies that Euclidean distance should be used.If
#'   \code{longlat} is \code{TRUE}, then the great circle
#'   distance is used to calculate the intercentroid
#'   distance.
#' @param type The type of scan statistic to implement.
#'   The default is \code{"poisson"}, with the
#'   other choice being \code{"binomial"}.
#' @param min.cases The minimum number of cases required for
#'   a cluster.  The default is 2.
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
#' @seealso \code{\link{scan.stat}},
#'   \code{\link{plot.scan}}, \code{\link{uls.test}},
#'   \code{\link{flex.test}}, \code{\link{dmst.test}},
#'   \code{\link{bn.test}}
#' @author Joshua French
#' @export
#' @references 
#' Kulldorff, M. (1997) A spatial
#'   scan statistic. Communications in Statistics -- Theory
#'   and Methods 26, 1481-1496.
#' 
#' Waller, L.A. and Gotway, C.A. (2005).
#'   Applied Spatial Statistics for Public Health Data.
#'   Hoboken, NJ: Wiley.  
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(longitude, latitude))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 alpha = 0.12, longlat = TRUE)
#' ## plot output for new york state
#' # specify desired argument values
#' mapargs = list(database = "state", region = "new york",
#' xlim = range(out$coords[,1]), ylim = range(out$coords[,2]))
#' # needed for "state" database (unless you execute library(maps))
#' data(stateMapEnv, package = "maps")
#' plot(out, usemap = TRUE, mapargs = mapargs)
#'
#' # a second example to match the results of Waller and Gotway (2005)
#' # in chapter 7 of their book (pp. 220-221).
#' # Note that the 'longitude' and 'latitude' used by them has
#' # been switched.  When giving their input to SatScan, the coords
#' # were given in the order 'longitude' and 'latitude'.
#' # However, the SatScan program takes coordinates in the order
#' # 'latitude' and 'longitude', so the results are slightly different
#' # from the example above.
#' coords = with(nydf, cbind(y, x))
#' out2 = scan.test(coords = coords, cases = floor(nydf$cases),
#'                   pop = nydf$pop, nsim = 0,
#'                   alpha = 1, longlat = TRUE)
#' # the cases observed for the clusters in Waller and Gotway: 117, 47, 44
#' # the second set of results match
#' c(out2$clusters[[1]]$cases, out2$clusters[[2]]$cases, out2$clusters[[3]]$cases)
scan.test = function(coords, cases, pop, 
                     ex = sum(cases)/sum(pop)*pop, 
                     nsim = 499, alpha = 0.1,  
                     ubpop = 0.5, longlat = FALSE, cl = NULL,
                     type = "poisson",
                     min.cases = 2) {
  # argument checking
  arg_check_scan_test(coords, cases, pop, ex, nsim, alpha, 
                      nsim + 1, ubpop, longlat, TRUE, 
                      k = 1, w = diag(nrow(coords)))
  if (length(min.cases) != 1 | min.cases < 1) {
    stop("min.cases must be a single number and >= 1")
  }
  
  # convert to proper format
  coords = as.matrix(coords)
  N = nrow(coords)
  # short names
  y = cases; e = ex
  # compute inter-centroid distances
  d = sp::spDists(coords, longlat = longlat)
  
  # for each region, determine sorted nearest neighbors
  # subject to population constraint
  nn = nnpop(d, pop, ubpop)

  # determine total number of cases in each successive 
  # window, total number of cases
  yin = nn.cumsum(nn, cases)
  ty = sum(y) # sum of all cases
  
  # compute test statistics for observed data
  if (type == "poisson") {
    ein = nn.cumsum(nn, ex)
    eout = ty - ein
    popin = NULL
    popout = NULL
    tpop = NULL
    tobs = stat.poisson(yin, ty - yin, ein, eout)
  } else if (type == "binomial") {
    ein = NULL
    eout = NULL
    tpop = sum(pop)
    popin = nn.cumsum(nn, pop)
    popout = tpop - popin
    tobs = stat.binom(yin, ty - yin, ty, popin, popout, tpop)
  }
  
  # determine distinct zones
  pri = randtoolbox::get.primes(N)
  wdup = duplicated(unlist(lapply(nn, function(x) cumsum(log(pri[x])))))

  # remove zones with a test statistic of 0 or don't have
  # min number of cases or are duplicted
  w0 = which(tobs == 0 | yin < min.cases | wdup)
  
  # determine zones
  zones = unlist(lapply(nn, function(x) sapply(seq_along(x), function(i) x[seq_len(i)])), recursive = FALSE)
  
  # remove zones with a test statistic of 0
  zones = zones[-w0]
  tobs = tobs[-w0]
  
  # compute test statistics for simulated data
  if (nsim > 1) {
    message("computing statistics for simulated data:")
    tsim = scan.sim(nsim = nsim, nn = nn, ty = ty,
                    ex = ex,
                    type = type, ein = ein, eout = eout,
                    popin = popin, popout = popout, tpop = tpop,
                    cl = cl)
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
  
  prep.scan(tobs = tobs, zones = zones, pvalue = pvalue, 
            coords = coords, cases = cases, pop = pop,
            ex = ex, longlat = longlat, w = NULL,
            d = d)
  # 
  # # order zones from largest to smallest test statistic
  # ozones = order(tobs, decreasing = TRUE)
  # zones = zones[ozones]
  # tobs = tobs[ozones]
  # 
  # # determine significant non-overlapping clusters
  # sig = smacpod::noc(zones)
  # 
  # # for the unique, non-overlapping clusters in order of significance,
  # # find the associated test statistic, p-value, centroid,
  # # window radius, cases in window, expected cases in window, 
  # # population in window, standarized mortality ration, relative risk,
  # sig_regions = zones[sig]
  # sig_tstat = tobs[sig]
  # sig_p = pvalue[ozones[sig]]
  # centroid = sapply(sig_regions, utils::head, n = 1)
  # boundary = sapply(sig_regions, utils::tail, n = 1)
  # sig_coords = coords[sapply(sig_regions, function(x) x[1]),, drop = FALSE]
  # sig_r = d[cbind(centroid, boundary)]
  # sig_yin = sapply(sig_regions, function(x) sum(y[x]))
  # sig_ein = sapply(sig_regions, function(x) sum(e[x]))
  # sig_popin = sapply(sig_regions, function(x) sum(pop[x]))
  # sig_smr = sig_yin/sig_ein
  # sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(sum(pop) - sig_popin))
  # sig_w = lapply(sig_regions, function(x) {
  #   matrix(c(0, rep(1, length(x) - 1)), nrow = 1)  
  # })
  # 
  # # reformat output for return
  # clusters = vector("list", length(sig_regions))
  # for (i in seq_along(clusters)) {
  #   clusters[[i]]$locids = sig_regions[[i]]
  #   clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
  #   clusters[[i]]$r = sig_r[i]
  #   clusters[[i]]$pop = sig_popin[i]
  #   clusters[[i]]$cases = sig_yin[i]
  #   clusters[[i]]$expected = sig_ein[i]
  #   clusters[[i]]$smr = sig_smr[i]
  #   clusters[[i]]$rr = sig_rr[i]
  #   clusters[[i]]$loglikrat = sig_tstat[[i]]
  #   clusters[[i]]$pvalue = sig_p[i]
  #   clusters[[i]]$w = sig_w[[i]]
  # }
  # outlist = list(clusters = clusters, coords = coords)
  # class(outlist) = "scan"
  # return(outlist)
}

# argument checking for all scan tests
arg_check_scan_test = 
  function(coords, cases, pop, ex, nsim, alpha, nreport,
           ubpop, longlat, parallel, k, w, type = NULL) {
  if(!(is.matrix(coords) | is.data.frame(coords))) stop("coords should be a matrix or a data frame")
  if(ncol(coords) != 2) stop("coords must have two columns")
  N = nrow(coords)
  if(length(cases) != N) stop("length(cases) != nrow(coords)")
  if(!is.numeric(cases)) stop("cases should be a numeric vector")
  if(length(pop) != N) stop("length(pop) != nrow(coords)")
  if(!is.numeric(pop)) stop("pop should be a numeric vector")
  if(length(ex) != N) stop("length(ex) != nrow(coords)")
  if(!is.numeric(ex)) stop("ex should be a numeric vector")
  if(length(alpha) != 1 || !is.numeric(alpha)) stop("alpha should be a numeric vector of length 1")
  if(alpha < 0 || alpha > 1) stop("alpha should be a value between 0 and 1")
  if(length(nsim) != 1 || !is.numeric(nsim)) stop("nsim should be a vector of length 1")
  if(nsim < 0) stop("nsim should be an non-negative integer")
  if(length(ubpop) != 1 || !is.numeric(ubpop)) stop("ubpop should be a numeric vector of length 1")
  if(ubpop<= 0 || ubpop > 1) stop("ubpop should be a value between 0 and 1")
  if(length(longlat) != 1) stop("length(longlat) != 1")
  if(!is.logical(longlat)) stop("longlat should be a logical value")
  if(length(parallel) != 1) stop("length(parallel) != 1")
  if(!is.logical(parallel)) stop("parallel should be a logical value")
  if(length(k) != 1) stop("k must have length 1")
  if(k < 1) stop("k must be an integer >= 1")
  if(!is.matrix(w)) stop("w must be a matrix")
  if(nrow(w) != ncol(w)) stop("w much be a square matrix")
  if(!is.numeric(w)) stop("w must be a numeric matrix")
  if(nrow(w) != nrow(coords)) stop("nrow(w) != nrow(coords)")
  if(floor(k) > nrow(coords)) stop("k cannot be more than the number of regions.")
  if (!is.null(type)) {
    if (!is.element(type, c("poisson", "binomial"))) {
      stop("type must be 'poisson' or 'binomial'")
    }
  }
}


