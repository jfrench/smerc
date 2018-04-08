#' Besag-Newell Test
#' 
#' \code{bn.test} implements the Besag-Newell test of 
#' Besag and Newell (1991) for finding disease clusters.
#' 
#' @param noc A logical value indicating whether all significant 
#' clusters should be returned (\code{FALSE}) or only the 
#' non-overlapping clusters (\code{TRUE}) arranged in 
#' order of significance.  The default is \code{TRUE}.
#' @param modified A logical value indicating whether a 
#' modified version of the test should be performed.  The 
#' original paper recommends computing the p-value for each
#' region as \code{1 - ppois(cstar - 1, lambda = expected)}.
#' The modified version replaces cstar with \code{yin}.
#' The default is \code{FALSE}.
#' @inheritParams scan.test
#' @inheritParams casewin
#'
#' @return Returns a list of length two of class \code{scan}. The first element (clusters) is a list containing the significant clusters and has the the following components: 
#' \item{locids}{The location ids of regions in a significant cluster.}
#' \item{coords}{The centroid of the initial region.}
#' \item{r}{The maximum radius of the cluster (in terms of intercentroid distance from the starting region).}
#' \item{pop}{The total population in the cluser window.}
#' \item{cases}{The observed number of cases in the cluster window.}
#' \item{expected}{The expected number of cases in the cluster window.}
#' \item{smr}{Standarized mortaility ratio (observed/expected) in the cluster window.}
#' \item{rr}{Relative risk in the cluster window.}
#' \item{tstat}{The loglikelihood ratio for the cluster window (i.e., the log of the test statistic).}
#' \item{pvalue}{The pvalue of the test statistic associated with the cluster window.}
#' \item{w}{The adjacency matrix of the cluster.}
#' The second element of the list is the centroid coordinates.  This is needed for plotting purposes.
#' @author Joshua French
#' @seealso \code{\link{scan.stat}}, \code{\link{plot.scan}}, 
#' \code{\link{scan.test}}, \code{\link{flex.test}}, 
#' \code{\link{dmst.test}}, \code{\link{uls.test}},
#' \code{\link{mlf.test}}
#' @export
#' @references Besag, J. and Newell, J.  (1991). The detection of clusters in rare diseases, Journal of the Royal Statistical Society, Series A, 154, 327-333.
#' @examples 
#' data(nydf)
#' data(nyw)
#' coords = with(nydf, cbind(x, y))
#' out = bn.test(coords = coords, cases = nydf$cases,
#'               pop = nydf$pop, cstar = 6, 
#'               alpha = 0.1)
#' plot(out)
#' 
#' data(nypoly)
#' library(sp)
#' plot(nypoly, col = color.clusters(out))
bn.test = function(coords, cases, pop, cstar,
                    alpha = 0.10,
                    lonlat = FALSE, noc = TRUE,
                    modified = FALSE) {
  # sanity checking
  arg_check_bn_test(coords, cases, pop, cstar, lonlat,
                    alpha, noc)

  coords = as.matrix(coords)
  N = nrow(coords)
  # estimate of constant risk
  r = sum(cases)/sum(pop)
  # intercentroid distances 
  d = sp::spDists(coords, longlat = lonlat)
  
  # find smallest windows with at least c* cases
  cwins = casewin(d, cases, cstar)
  # determine size of each window
  l = sapply(cwins, length)
  # determine population in each window
  case_cwins = sapply(cwins, function(x) sum(cases[x]))
  ncwins = sapply(cwins, function(x) sum(pop[x]))
  ex = r * ncwins
  
  if (!modified) {
    # pvalue = stats::ppois(cstar - 1, lambda = ex, lower.tail = FALSE)
    pvalue = stats::ppois(cstar, lambda = ex, lower.tail = FALSE)
    } else {
    pvalue = stats::ppois(case_cwins - 1, lambda = ex, lower.tail = FALSE)
  }

  op = order(pvalue)

  if (noc) {
    # determine idx of unique non-overlapping clusters in
    # order of significance
    u = smacpod::noc(cwins[op])
    op = op[u]
    # return only significant clusters
    if(pvalue[op][1] > alpha) {
      warning("No significant clusters.  Returning most likely cluster.")
      op = op[1]
    } else {
      op = op[which(pvalue[op] <= alpha)]
    }
  }
 
  # for the unique, non-overlapping clusters in order of significance,
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ratio, 
  # relative risk,
  sig_regions = cwins[op]
  sig_tstat = l[op]
  sig_p = pvalue[op]
  sig_coords = coords[op,, drop = FALSE]
  sig_r = d[op, l[op]]
  sig_yin = case_cwins[op]
  sig_ein = ex[op]
  sig_popin = ncwins[op]
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((sum(cases) - sig_yin)/(sum(pop) - sig_popin))
  sig_w = lapply(sig_tstat, function(x) matrix(c(0, rep(1, x - 1)), nrow = 1))

  # reformat output for return
  clusters = vector("list", length(op))
  for (i in seq_along(clusters))   {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
    clusters[[i]]$r = sig_r[i]
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$tstat = sig_tstat[i]
    clusters[[i]]$pvalue = sig_p[i]
    clusters[[i]]$w = sig_w[[i]]
  }
  outlist = list(clusters = clusters, coords = coords)
  class(outlist) = "scan"
  return(outlist)
}

arg_check_bn_test = 
  function(coords, cases, pop, cstar, lonlat, alpha, noc)
  {
    if(!(is.matrix(coords) | is.data.frame(coords))) stop("coords should be a matrix or a data frame")
    if(ncol(coords) != 2) stop("coords must have two columns")
    N = nrow(coords)
    if(length(cases) != N) stop("length(cases) != nrow(coords)")
    if(!is.numeric(cases)) stop("cases should be a numeric vector")
    if(length(pop) != N) stop("length(pop) != nrow(coords)")
    if(!is.numeric(pop)) stop("pop should be a numeric vector")
    if(length(cstar) != 1 || !is.numeric(cstar)) stop("cstar should be a numeric vector of length 1")
    if(cstar < 1 || cstar > sum(cases)) stop("cstar should be at least 1 and less than or equal to the sum(cases)")
    if(length(lonlat) != 1) stop("length(lonlat) != 1")
    if(!is.logical(lonlat)) stop("lonlat should be a logical value")
    if(length(alpha) != 1 || !is.numeric(alpha)) stop("alpha should be a numeric vector of length 1")
    if(alpha < 0 || alpha > 1) stop("alpha should be a value between 0 and 1")
    if(length(noc) != 1) stop("length(noc) != 1")
    if(!is.logical(noc)) stop("noc should be a logical value")
}
