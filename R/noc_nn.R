#' Returned ordered non-overlapping clusters
#'
#' \code{noc_nn} returns the non-overlapping clusters
#' ordered by the values on \code{tobs_nn}.
#' @param nn A nearest neighbor list.
#' @param tobs_nn The test statistics for each sequence
#'   candidate zone.
#' @param nnoc The number of non-overlapping clusters to
#'   return.
#'
#' @return A list with \code{clusts}, a list with elements
#'   containing the non-overlapping clusters, and
#'   \code{tobs}, the associated test statistic for each
#'   cluster.
#' @export
#' @keywords internal
noc_nn <- function(nn, tobs_nn, nnoc = length(nn)) {
  # fix bug when ubpop is smaller than proportion of maximum
  # population region
  remaining <- which(lapply(nn, length) > 0) # seq_along(nn)
  clusts <- list()
  tobs_noc <- list()
  current_tobs_max <- 1
  starting_region_list <- unlist(lapply(nn, `[`, i = 1))

  while (length(remaining) > 0 &
    length(clusts) < nnoc &
    current_tobs_max > 0) {
    # find the max of tobs_nn and corresponding index for
    # each nn row
    tobs_nn_max_idx <- unlist(lapply(tobs_nn[remaining], which.max))
    tobs_nn_max <- unlist(lapply(tobs_nn[remaining], max))

    # find the max over tobs_max_nn and corresponding index
    # mlc <- max(unlist(tobs_nn_max))
    mlc_max_idx <- which.max(tobs_nn_max)
    mlc_max_seq_idx <- which.max(tobs_nn[remaining][[mlc_max_idx]])
    clusts <- append(clusts, list(nn[remaining][[mlc_max_idx]][seq_len(mlc_max_seq_idx)]))
    current_tobs_max <- tobs_nn_max[mlc_max_idx]
    tobs_noc <- append(tobs_noc, list(current_tobs_max))
    mlc_max_seq_idx <- which.max(tobs_nn[remaining][[mlc_max_idx]])
    clusts_unlist <- unlist(clusts)
    remaining <- setdiff(remaining, clusts_unlist)
    # remaining <- which(!is.element(starting_region_list, clusts_unlist))
    # tobs_nn[clusts_unlist] <- numeric(0)
    # prune tobs_nn
    for (i in remaining) {
      wmin <- which(is.element(nn[[i]], clusts_unlist))
      # if there's at least one overlapping regions,
      # prune the observed statistics
      if (length(wmin) > 0) {
        tobs_nn[[i]] <- tobs_nn[[i]][seq_len(wmin[1] - 1)]
      }
    }
  }
  return(list(
    clusts = clusts,
    tobs = unlist(tobs_noc, use.names = FALSE)
  ))
}
