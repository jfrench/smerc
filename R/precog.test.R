#' PreCoG Scan Test
#'
#' \code{precog.test} is an implementation of the
#' Prefiltered Component-based Greedy Scan Method.
#' @inheritParams scan.test
#' @inheritParams rflex.test
#' @param tol_prob A single numeric value between 0 and 1
#'   that describes the quantile of the tolerance envelopes
#'   used to prefilter regions from the candidate zones.
#' @param ysim A matrix of size \code{nsim}\eqn{\times n},
#'   where \eqn{n} is the number of regions in the study
#'   area. This is a matrix of \code{nsim} realizations of
#'   the case counts for each region in the study area under
#'   the null hypothesis. This argument is only not meant to
#'   be used by the user.
#'
#' @return Returns a \code{smerc_cluster} object.
#' @seealso \code{\link{print.smerc_cluster}},
#'   \code{\link{summary.smerc_cluster}},
#'   \code{\link{plot.smerc_cluster}},
#' @author Joshua French
#' @export
#' @examples
#' data(nydf)
#' data(nyw)
#' out <- precog.test(coords = nydf[,c("x", "y")],
#'                    cases = floor(nydf$cases),
#'                    pop = nydf$pop, w = nyw, nsim = 19,
#'                    alpha = 0.2)
#' plot(out)
precog.test <- function(coords,
                        cases,
                        pop,
                        w,
                        ex = sum(cases) / sum(pop) * pop,
                        nsim = 499,
                        tol_prob = 0.90,
                        alpha = 0.1,
                        ubpop = 0.5,
                        longlat = FALSE,
                        cl = NULL,
                        ysim = NULL) {

  arg_check_scan_test(
    coords = coords, cases = cases,
    pop = pop, ex = ex, nsim = nsim,
    alpha = alpha, ubpop = ubpop,
    longlat = longlat,
    k = 1, w = w,
    type = "poisson", simdist = "multinomial",
    min.cases = 1
  )
  arg_check_tol_prob(tol_prob)
  if (nsim < 1) {
    stop("nsim must be at least 1 (and should probably be much larger)")
  }

  # compute preliminary data
  ty <- sum(cases) # total number of cases
  max_pop <- ubpop * sum(pop) # max population allowed
  # distance matrix
  d <- gedist(as.matrix(coords), longlat = longlat)
  # initial zones
  zones <- apply((w + diag(nrow(w))) == 1, 1, which,
                 simplify = FALSE)
  # expected in each candidate
  ein  <-  zones.sum(zones, ex)
  logein <- log(ein)
  eout <- ty - ein
  logeout <- log(eout)

  # determine cases in initial candidate zones
  yin <- zones.sum(zones, cases)
  # associated test statistics
  tobs <- stat_poisson_adj(yin, ty, logein, logeout)

  # determine tolerance envelopes for null data based on
  # connected neighbor candidate zones
  message("computing statistics for simulated data:")
  precog_null_info <-
    precog.sim(nsim = nsim, zones = zones, ty = ty, ex = ex,
               w = w, pop = pop, max_pop = max_pop,
               logein = logein, logeout = logeout, d = d,
               cl = cl, tol_prob = tol_prob, ysim = ysim)

  # keep regions with test statistic above tol_enveloeps
  # keep <- unlist(cc, use.names = FALSE)
  keep <- which(tobs > precog_null_info$tol_envelopes)
  # w for connected clusters
  ccw <- w[keep, keep]

  # for each connected cluster, run dmst method
  tall <- mst.all(nndist(d[keep, keep], 1),
                  cases = cases[keep], pop = pop[keep],
                  w = ccw,
                  ex = ex[keep], ty = ty,
                  max_pop = max_pop, type = "all",
                  early = FALSE, nlinks = "one",
                  progress = FALSE)

  # extract sequential nearest neighbors for each starting region
  nn2 <- lgetElement(tall, name = "locids")
  # extract associated test statistics of nn2
  tobs_nn <- lgetElement(tall, name = "loglikrat")

  # order candidate zones by largest test statistic and
  # and prune overlapping clusters
  noc_info <- noc_nn(nn2, tobs_nn)

  # convert cc cluster region ids to original region ids
  noc_info$clusts <- lapply(noc_info$clusts, function(x) keep[x])
  # test statistic associated with each cluster
  tobs <- noc_info$tobs

  # compute pvalue
  tsim <- precog_null_info[[1]]
  pvalue <- mc_pvalue_cpp(tobs, tsim)

  # prune insignificant clusters
  pruned <- sig_prune(tobs = tobs,
                      zones = noc_info$clusts,
                      pvalue = pvalue,
                      alpha = alpha)
  smerc_cluster(
    tobs = pruned$tobs, zones = pruned$zones,
    pvalue = pruned$pvalue, coords = coords,
    cases = cases, pop = pop, ex = ex,
    longlat = longlat, method = "precog scan",
    rel_param = list(
      type = "poisson",
      simdist = "multinomial",
      nsim = nsim,
      ubpop = ubpop,
      tol_prob = tol_prob
    ),
    alpha = alpha,
    w = w, d = d
  )
}
