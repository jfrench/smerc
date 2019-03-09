#' Construct connected subgraphs
#'
#' \code{csg}, \code{lcsg}, and \code{scsg} construct
#' connected subgraphs.
#' \code{set} contains a vector of vertices of the pattern 1, 2, 3, ..., N.
#' \code{idx} is a vector of possible vertices being considered as a subgraph.
#' \code{w} is a connectivity matrix relating the N vertices.
#' \code{w[i,j] = 1} if vertices i and j are connected, i.e., if they share an edge.
#' The dimensions of \code{w} are \eqn{N times k}, where \code{k = length(idx)}.
#' While the rows of \code{w} contain adjacency information for all N
#' vertices, only the \code{idx} columns of the complete adjacency matrix
#' are used in \code{w}.  See Details for discussion of \code{scsg}.
#'
#' \code{scsg} performs
#' a sequence of \code{lcsg} calls.  Starting with \code{lset == list(idx[1])},
#' \code{scsg} keeps iteratively building more connected subsgraphs by perfoming
#' something like:  set1 = list(idx[1]).  set2 = lcsg(set1, idx, w).
#' set3 = lcsg(set2, idx, w).  This is done until there are no more connected
#' subgraphs among the elements of \code{idx}.
#'
#' @param set A vector of (presumably connected) vertices.
#' @param lset A list of sets.
#' @param idx A vector of vertices considered for inclusion in the
#' subgraph, e.g., based on nearest neighbors.
#' @param w The adjacency matrix for all vertices by row, but with only the \code{idx} columns
#' @param verbose A logical value indicating whether very descriptive messages
#' should be provided.  Default is \code{FALSE}.  If \code{TRUE}, this can
#' be useful for diagnosing where the sequences of connected subgraphs is
#' slowing down/having problems.
#'
#' @return A list of with all possible combinations of \code{set} and
#' each possible connected vertex in \code{idx}, or \code{NULL} if none
#' are possible.
#' @export
#'
#' @examples
#' data(nydf)
#' data(nyw)
#' # determine 50 nn of region 1 for NY data
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' d = sp::spDists(coords, longlat = TRUE)
#' nn50 = order(d[1,])[1:50]
#' w = nyw[,nn50]
#' set = 1
#' # first set of connected neighbors
#' nb1 = csg(set, idx = nn50, w = w)
#' # extend set of connected neighbors
#' # for first element of nb1
#' set2 = nb1[[1]]
#' nb2 = csg(set2, idx = nn50, w = w)
#' # do the same thing for all sets in nb1
#' nb2e = lcsg(nb1, idx = nn50, w = w)
#' # the sets in nb2 should be present in the
#' # first 9 positions of nb2e
#' all.equal(nb2, nb2e[seq_along(nb2)])
#'
#' # apply scsg to first 10 nn of vertex 1
#' nn10 = order(d[1,])[1:10]
#' w = nyw[, nn10]
#' nb3 = scsg(nn10, w, verbose = TRUE)
csg = function(set, idx, w) {
  if (ncol(w) != length(idx)) {
    stop("ncol(w) != length(idx)")
  }
  # determine neighbors of current set
  nb = setdiff(idx[which(matrixStats::colMaxs(w[set,, drop = FALSE]) == 1)],
               set)
  # if there are some neighbors
  if (length(nb) > 0) {
    # add neighbors (individually) to current set
    newsets = lapply(nb, function(x) c(set, x))
    # return unique sets
    return(newsets[distinct(newsets)])
  } else {
    # return NULL if there are no neighbors
    return(NULL)
  }
}

#' @export
#' @rdname csg
lcsg = function(lset, idx, w) {
  x = unlist(lapply(lset, csg, w = w, idx = idx),
             recursive = FALSE)
  if (is.null(x)) {
    return(NULL)
  } else {
    x[distinct(x)]
  }
}

#' @export
#' @rdname csg
scsg = function(idx, w, verbose = FALSE) {
  # initial set
  nidx = length(idx)
  out = vector("list", nidx)
  out[[1]] = list(idx[1])

  # set stopping conditions
  # stop if nidx == 1, otherwise, consider expanding subgraph
  stop = (nidx == 1)
  j = 2
  if (verbose) {
    txt = paste("1/", nidx, ". Adding region ", idx[1],
                " at ", Sys.time(), ".", sep = "")
    message(txt)
  }
  while (!stop) {
    if (verbose) {
      txt = paste(j, "/", nidx, ". Adding region ", idx[j],
                  " at ", Sys.time(), ".", sep = "")
      message(txt)
    }
    # successively generate list of connected subgraphs
    temp = lcsg(out[[j - 1]], idx = idx, w = w)
    # stop the loop of there are no more connected subgraphs
    # to add (i.e., out[[j]] is NULL)
    if (is.null(temp) | j == nidx) {
      stop = TRUE
    }
    out[[j]] = temp
    j = j + 1
  }
  return(unlist(out, recursive = FALSE))
}
