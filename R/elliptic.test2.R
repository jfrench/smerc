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
#'
#' @inherit scan.test return params
#' @seealso \code{\link{print.smerc_cluster}},
#' \code{\link{summary.smerc_cluster}},
#' \code{\link{plot.smerc_cluster}},
#' \code{\link{scan.stat}}, \code{\link{scan.test}}
#' @author Joshua French
#' @export
#' @references Kulldorff, M. (1997) A spatial scan
#'   statistic. Communications in Statistics - Theory and
#'   Methods, 26(6): 1481-1496,
#'   <doi:10.1080/03610929708831995>
#'
#' Kulldorff, M., Huang, L., Pickle,
#' L. and Duczmal, L. (2006) An elliptic spatial scan
#' statistic. Statististics in Medicine, 25:3929-3943.
#' <doi:10.1002/sim.2490>
#' @examples
#' data(nydf)
#' coords <- nydf[, c("x", "y")]
#' out <- elliptic.test(
#'   coords = coords,
#'   cases = floor(nydf$cases),
#'   pop = nydf$pop, ubpop = 0.1,
#'   nsim = 2,
#'   alpha = 0.12
#' )
elliptic.test <- function(coords, cases, pop,
                          ex = sum(cases) / sum(pop) * pop,
                          nsim = 499, alpha = 0.1,
                          ubpop = 0.5,
                          shape = c(1, 1.5, 2, 3, 4, 5),
                          nangle = c(1, 4, 6, 9, 12, 15),
                          a = 0.5,
                          cl = NULL,
                          type = "poisson",
                          min.cases = 2) {
  # argument checking
  arg_check_scan_test(
    coords = coords, cases = cases,
    pop = pop, ex = ex, nsim = nsim,
    alpha = alpha,
    ubpop = ubpop, longlat = FALSE,
    parallel = FALSE, k = 1,
    w = diag(nrow(coords)),
    type = type, min.cases = min.cases
  )
  arg_check_elliptic_test(shape = shape, nangle = nangle, a = a)

  # convert to proper format
  coords <- as.matrix(coords)
  N <- nrow(coords)

  enn <- elliptic.nn(coords,
    pop = pop, ubpop = ubpop,
    shape = shape, nangle = nangle
  )
  # determine the expected cases in each successive
  # zone, on log scale, total number of cases,
  # log of expected cases outside zone
  ein <- nn.cumsum(enn$nn, ex)
  logein <- log(ein)
  ty <- sum(cases) # sum of all cases
  logeout <- log(ty - ein)

  # determine yin and yout for all windows for observed data
  yin <- nn.cumsum(enn$nn, cases)

  ### calculate scan statistics for observed data
  # of distance from observation centroid
  pen <- elliptic.penalty(a, enn$shape_all)
  tobs <- stat_poisson_adj(
    yin = yin, ty = ty,
    logein = logein, logeout = logeout,
    a = a, pen = pen,
    min.cases = min.cases
  )

  # tobs in nn format
  nnn <- unlist(lapply(enn$nn, length), use.names = FALSE)
  tobs_nn <- split(tobs, f = rep(seq_along(enn$nn), times = nnn))

  # determine non-overlapping clusters in order of significance
  noc_info <- noc_enn(enn$nn, tobs_nn, enn$shape_nn, enn$angle_nn)
  tobs <- noc_info$tobs

  # setup list for call
  if (nsim > 0) {
    tsim <- elliptic.sim.adj(
      nsim = nsim, ex = ex,
      nn = enn$nn, ty = ty,
      logein = logein, logeout = logeout,
      a = a, pen = pen, min.cases = min.cases,
      cl = cl
    )

    # p-values associated with these max statistics for each centroid
    pvalue <- mc.pvalue(tobs, tsim)
  } else {
    pvalue <- rep(1, length(tobs))
  }

  # # construct all zones
  # zones = nn2zones(enn$nn)
  # zones = zones[-w0]

  # significant, ordered, non-overlapping clusters and
  # information
  # pruned = sig_noc(tobs = tobs, zones = zones,
  #                  pvalue = pvalue, alpha = alpha,
  #                  order_by = "tobs")
  # significant, ordered, non-overlapping clusters and
  # information
  pruned <- sig_prune(
    tobs = tobs, zones = noc_info$clusts,
    pvalue = pvalue, alpha = alpha
  )
  shape_all <- noc_info$shape[seq_along(pruned$tobs)]
  angle_all <- noc_info$angle[seq_along(pruned$tobs)]

  smerc_cluster(
    tobs = pruned$tobs, zones = pruned$zones,
    pvalue = pruned$pvalue, coords = coords,
    cases = cases, pop = pop, ex = ex,
    longlat = FALSE, method = "elliptic",
    rel_param = list(
      type = type,
      simdist = "multinomial",
      nsim = nsim,
      ubpop = ubpop,
      min.cases = min.cases,
      a_penalty = a,
      shapes = shape,
      nangles = nangle
    ),
    alpha = alpha,
    w = NULL, d = NULL, a = a,
    shape_all = shape_all,
    angle_all = angle_all
  )
}

#' Additional argument checking for elliptic_test
#'
#' @param shape A vector of shapes (values >= 1)
#' @param nangle A vector of angles for each shape (values >= 1)
#' @param a A penalty parameter (a >= 0)
#' @return NULL
#' @noRd
arg_check_elliptic_test <-
  function(shape, nangle, a) {
    if (length(shape) != length(nangle)) {
      stop("The length of shape and nangle must match.")
    }
    arg_check_shape(shape)
    arg_check_nangle(nangle)
    arg_check_a(a)
  }
