#' Additional argument checking for elliptic_test
#'
#' @param shape A vector of shapes (values >= 1)
#' @param nangle A vector of angles for each shape (values >= 1)
#' @param a A penalty parameter (a >= 0)
#' @return NULL
#' @noRd
arg_check_elliptic_test =
  function(shape, nangle, a) {
    if (length(shape) != length(nangle)) {
      stop("The length of shape and nangle must match.")
    }
    arg_check_shape(shape)
    arg_check_nangle(nangle)
    arg_check_a(a)
  }
