#' @method clusters smerc_cluster
#' @export
clusters.smerc_cluster = function(x, idx = seq_along(x$clusters), ...) {
  if (min(idx) < 1 | max(idx) > length(x$clusters)) {
    stop("invalid idx values")
  }
  lget(x$clusters[idx], "locids")
}
