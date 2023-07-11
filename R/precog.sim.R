precog.sim <- function(nsim = 1, zones, ty, ex, w, pop,
                       max_pop, logein, logeout, d,
                       cl = NULL, tol_prob = 0.90,
                       ysim = NULL) {

  if (is.null(ysim)) {
    ysim <- t(stats::rmultinom(nsim, size = ty, prob = ex))
  }

  # get statistics for each null data set
  stats_null <-
    apply(
      ysim,
      MARGIN = 1,
      FUN = function(simcases) {
        yin <- zones.sum(zones, simcases)
        stat_poisson_adj(yin, ty, logein, logeout, min.cases = 1)
      })

  # get tolerance envelopes for each region
  tol_envelopes <- apply(stats_null, 1, stats::quantile,
                         prob = tol_prob)

  # filter regions for null data sets
  cc_keep <- apply(stats_null, 2, function(tobs) {
    which(tobs > tol_envelopes)
  }, simplify = FALSE)

  # determine maximum test statistic for each simulated
  # data set
  max_stats2_null <- pbapply::pblapply(seq_along(cc_keep),
                                       function(i) {
    keep <- cc_keep[[i]]
    tall <- mst.all(nndist(d[keep, keep], 1),
                    cases = ysim[i, keep], pop = pop[keep],
                    w = w[keep, keep],
                    ex = ex[keep], ty = ty,
                    max_pop = max_pop, type = "maxonly",
                    early = FALSE, nlinks = "one",
                    progress = FALSE
    )
    return(max(tall))
  }, cl = cl)
  # reformat test statistics
  tsim <- unlist(max_stats2_null, use.names = FALSE)
  # return simulation information
  return(list(tsim = tsim,
              tol_envelopes = tol_envelopes))
}