#' Argument checking for *.sim functions
#'
#' Check the arguments of the \code{*.sim} functions.
#' @return NULL
#' @export
#' @keywords internal
arg_check_sim = function(nsim,
                         ty,
                         ex,
                         type,
                         nn = NULL,
                         zones = NULL,
                         ein = NULL,
                         eout = NULL,
                         tpop = NULL,
                         popin = NULL,
                         popout = NULL,
                         w = NULL,
                         pop = NULL,
                         ubpop = NULL,
                         static) {
  if (length(nsim) != 1 | !is.numeric(nsim) | nsim < 1) {
    stop("nsim must be a positive integer")
  }
  if (length(ty) != 1) stop("ty must be a single number")
  if (!is.numeric(ex)) stop("ex must be a numeric vector")
  if (length(type) != 1 |
      !is.element(type, c("poisson", "binomial"))) {
    stop("type must be 'poisson' or 'binomial'")
  }
  if (static) {
    if (!is.null(nn)) {
      if (!is.list(nn)) stop("nn must be a list")
      nz = sum(sapply(nn, length))
    }
    if (!is.null(zones)) {
      if (!is.list(zones)) stop("zones must be a list")
      nz = length(zones)
    }
    if (type == "poisson") {
      if (is.null(ein) | is.null(eout)) {
        stop("ein and eout must be provided when type is 'poisson'")
      }
      if (nz != length(ein)) {
        stop("ein has improper length")
      }
      if (nz != length(eout)) {
        stop("eout has improper length")
      }
    }
    if (type == "binomial") {
      if (is.null(popin) | is.null(popout) | is.null(tpop)) {
        stop("popin, popout, and tpop must be provided when type is 'binomial'")
      }
      if (nz != length(popin)) {
        stop("popin has improper length")
      }
      if (nz != length(popout)) {
        stop("popout has improper length")
      }
      if (length(tpop) != 1) {
        stop("tpop must be a single number")
      }
    }
  }
  if (!is.null(w)) {
    if (is.null(dim(w))) stop("w must be matrix-like")
    if (length(dim(w)) != 2) stop("w must be two-dimensional")
    if (nrow(w) != ncol(w)) stop("w must be square")
  }
  if (!is.null(tpop)) {
    if (length(tpop) != 1 | !is.numeric(tpop)) {
      stop("tpop must be a single number")
    }
  }
  if (!is.null(w)) {
    if (is.null(dim(w))) stop("w must be matrix-like")
    if (length(dim(w)) != 2) stop("w must be two-dimensional")
    if (nrow(w) != ncol(w)) stop("w must be square")
  }
  if (!is.null(pop)) {
    if (!is.numeric(pop) | length(pop) != length(ex)) {
      stop("pop must be a numeric vector with length = length(ex)")
    }
  }
  if (!is.null(ubpop)) {
    if (length(ubpop) != 1 | ubpop <= 0 | ubpop > 1) {
      stop("ubpop must be in (0, 1]")
    }
  }
}
