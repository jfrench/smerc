#' Prune significant, non-overlapping zones
#'
#' \code{sig_prune} returns the significant
#' zones ordered from most significant to least significant (assuming the zones are already in order)
#' @inheritParams smerc_cluster
#'
#' @return A list with the significant, ordered,
#' non-overlapping \code{tobs}, \code{zones}, \code{pvalue}.,
#' and \code{idx} (a vector with the relevant indices of
#' the original zones).
#' @export
#' @keywords internal
#' @examples
#' tobs <- c(3, 2, 1)
#' zones <- list(1:2, 3:4, 5:6)
#' pvalue <- c(0.001, 0.05, 0.15)
#' sig_prune(tobs, zones, pvalue, alpha = 0.05)
sig_prune <- function(tobs, zones, pvalue, alpha) {
  # argument checking
  N <- length(tobs)
  arg_check_tobs(tobs)
  arg_check_zones(zones, N)
  arg_check_pvalue(pvalue, N)
  arg_check_alpha(alpha)

  # determine if there are any significant zones
  minp <- which.min(pvalue)
  if (pvalue[minp] > alpha) {
    warning("No significant clusters.  Returning most likely cluster.")
    sig <- minp
  } else {
    sig <- which(pvalue <= alpha)
  }
  # only keep significant zones and info
  return(list(
    tobs = tobs[sig],
    zones = zones[sig],
    pvalue = pvalue[sig]
  ))
}
