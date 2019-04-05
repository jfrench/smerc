#' Distance-based weights
#'
#' \code{dweights} constructs a distance-based weights matrix.
#' The \code{dweights} function can be used to
#' construct a weights matrix \code{w} using the method
#' of Tango (1995), Rogerson (1999), or a basic style.
#'
#' \code{coords} is used to construct an \eqn{n \times n}
#' distance matrix \code{d}.
#'
#' If \code{type = "basic"}, then
#' \eqn{w_{ij} = exp(-d_{ij}/\kappa)}.
#'
#' If \code{type = "rogerson"}, then
#' \eqn{w_{ij} = exp(-d_{ij}/\kappa)/\sqrt(pop_i/pop * pop_j/pop)}.
#'
#' If \code{type = "tango"}, then
#' \eqn{w_{ij} = exp(-4 * d_{ij}^2/\kappa^2)}.
#'
#'
#' @inheritParams scan.test
#' @param kappa A positive constant related to strength of
#' spatial autocorrelation.
#' @param type The type of weights matrix to construct.
#' Current options are \code{"basic"}, \code{"tango"},
#' and \code{"rogerson"}.  Default is \code{"basic"}.  See Details.
#'
#' @return Returns an \eqn{n \times n} matrix of weights.
#' @references
#' Tango, T.  (1995) A class of tests for detecting "general" and "focused" clustering of rare diseases.  Statistics in Medicine.  14:2323-2334.
#'
#' Rogerson, P. (1999) The Detection of Clusters Using A Spatial Version of the Chi-Square
#' Goodness-of-fit Test.  Geographical Analysis. 31:130-147
#' @author Joshua French
#' @seealso \code{\link{tango.test}}
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("longitude", "latitude")])
#' w = dweights(coords, kappa = 1, longlat = TRUE)
dweights = function(coords, kappa = 1, longlat = FALSE, type = "basic", pop = NULL) {
  arg_check_dweights(coords, kappa, longlat, type, pop)
  d = sp::spDists(as.matrix(coords), longlat = longlat)
  if (type == "basic") {
    w = exp(-d / kappa)
  } else if (type == "rogerson") {
    popp = pop / sum(pop)
    w = exp(-d / kappa) / sqrt(outer(popp, popp))
  } else if (type == "tango") {
    w = exp(-4 * (d / kappa) ^ 2)
  }
  return(w)
}

arg_check_dweights = function(coords, kappa, longlat, type, pop) {
  if (!(is.matrix(coords) | is.data.frame(coords))) {
    stop("coords should be a matrix or a data frame")
  }
  if (ncol(coords) != 2) {
    stop("coords must have two columns")
  }
  N = nrow(coords)
  if (length(kappa) != 1 || !is.numeric(kappa)) {
    stop("kappa should be a numeric vector of length 1")
  }
  if (kappa <= 0) stop("kappa must be positive")
  if (length(longlat) != 1) stop("length(longlat) != 1")
  if (!is.logical(longlat)) {
    stop("longlat should be a logical value")
  }
  if (length(type) != 1) {
    stop("type must be a single character")
  }
  if (!is.element(type, c("basic", "rogerson", "tango"))) {
    stop("invalid type")
  }
  if (!is.null(pop)) {
    if (length(pop) != N) {
      stop("length(pop) != nrow(coords)")
    }
    if (!is.numeric(pop)) {
      stop("pop should be a numeric vector")
    }
  }
  if (type == "rogerson") {
    if (is.null(pop)) {
      stop("pop must be provided when type = 'rogerson'")
    }
    if (any(pop == 0)) {
      stop("regions cannot contain zero population when type = 'rogerson'")
    }
  }
}
