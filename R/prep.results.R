#' Prepare scan results
#' 
#' \code{prep.scan} prepares the results of a scan-related
#' test and turns them into a scan object.
#'
#' @param tobs The observed test statistics for each zone
#' @param zones A list of zones
#' @param pvalue The p-value associated with each test statistic
#' @inheritParams flex.test
#' @param d A precomputed distance matrix based on coords
#'
#' @keywords internal
#' @return A scan object
#' @export
#'
#' @examples
prep.scan = function(tobs, zones, pvalue, coords, cases, 
                     pop, ex, longlat, w = NULL, d = NULL) {
  # order zones from largest to smallest test statistic
  ozones = order(tobs, decreasing = TRUE)
  zones = zones[ozones]
  tobs = tobs[ozones]
  
  # determine significant non-overlapping clusters
  sig = smacpod::noc(zones)
  
  # total cases and populatoin
  ty = sum(cases)
  tpop = sum(pop)
  
  # for the unique, non-overlapping clusters in order of significance,
  # find the associated test statistic, p-value, centroid,
  # window radius, cases in window, expected cases in window, 
  # population in window, standarized mortality ration, relative risk,
  sig_regions = zones[sig]
  sig_tstat = tobs[sig]
  sig_p = pvalue[ozones[sig]]
  centroid = sapply(sig_regions, utils::head, n = 1)
  boundary = sapply(sig_regions, utils::tail, n = 1)
  if (!is.null(d)) {
    sig_r = d[cbind(centroid, boundary)]
  } else {
    sig_r = rep(NA, length(centroid))
  }
  sig_max_dist = unname(sapply(sig_regions, function(x) {
    max(sp::spDists(coords[x, ,drop = FALSE]), longlat = longlat)
  }))
  sig_coords = coords[sapply(sig_regions, function(x) x[1]),, drop = FALSE]
  sig_yin = sapply(sig_regions, function(x) sum(cases[x]))
  sig_ein = sapply(sig_regions, function(x) sum(ex[x]))
  sig_popin = sapply(sig_regions, function(x) sum(pop[x]))
  sig_smr = sig_yin/sig_ein
  sig_rr = (sig_yin/sig_popin)/((ty - sig_yin)/(tpop - sig_popin))
  if (is.null(w)) {
    sig_w = lapply(sig_regions, function(x) {
      matrix(c(0, rep(1, length(x) - 1)), nrow = 1)  
    })
  } else {
    sig_w = sapply(sig_regions, function(x) w[x, x, drop = FALSE], simplify = FALSE)  
  }

  # reformat output for return
  clusters = vector("list", length(sig_regions))
  for (i in seq_along(clusters)) {
    clusters[[i]]$locids = sig_regions[[i]]
    clusters[[i]]$coords = sig_coords[i,, drop = FALSE]
    clusters[[i]]$r = sig_r[i]
    clusters[[i]]$max_dist = sig_max_dist[i]
    clusters[[i]]$pop = sig_popin[i]
    clusters[[i]]$cases = sig_yin[i]
    clusters[[i]]$expected = sig_ein[i]
    clusters[[i]]$smr = sig_smr[i]
    clusters[[i]]$rr = sig_rr[i]
    clusters[[i]]$loglikrat = sig_tstat[[i]]
    clusters[[i]]$pvalue = sig_p[i]
    clusters[[i]]$w = sig_w[[i]]
  }
  outlist = list(clusters = clusters, 
                 coords = coords,
                 number_of_locations = length(cases),
                 total_population = tpop,
                 total_cases = ty,
                 cases_per_100k = ty/tpop * 1e5)
  class(outlist) = "scan"
  return(outlist)
}