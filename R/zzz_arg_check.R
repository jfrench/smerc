#' Check coords argument
#'
#' @param coords An matrix-like object with N rows
#' @noRd
arg_check_coords = function(coords) {
  if (!(is.matrix(coords) | is.data.frame(coords))) {
    stop("coords should be a matrix or a data frame")
  }
  if (ncol(coords) != 2) {
    stop("coords must have two columns")
  }
}

#' Check cases argumnent
#'
#' @param cases A numeric vector of cases
#' @param N The number of rows in coords
#' @noRd
arg_check_cases = function(cases, N) {
  if (length(cases) != N) {
    stop("length(cases) != nrow(coords)")
  }
  if (!is.numeric(cases)) {
    stop("cases should be a numeric values")
  }
  if (!is.vector(cases)) {
    stop("cases should be a vector")
  }
  if (min(cases) < 0) {
    stop("cases must have non-negative values")
  }
}

#' Check population argument
#'
#' @param pop A vector of population values
#' @param N The number of rows in coords
#' @noRd
arg_check_pop = function(pop, N) {
  if (length(pop) != N) {
    stop("length(pop) != nrow(coords)")
  }
  if (!is.numeric(pop)) {
    stop("pop should be numeric values")
  }
  if (!is.vector(pop)) {
    stop("pop should be a vector")
  }
  if (min(pop) < 1) {
    stop("pop values must be >= 1")
  }
}

#' Check cstar argument
#'
#' @param cstar Case radius
#' @param cases Numeric vector of cases
#' @noRd
arg_check_cstar = function(cstar, cases) {
  if (length(cstar) != 1 || !is.numeric(cstar)) {
    stop("cstar should be a numeric vector of length 1")
  }
  if (cstar < 1 || cstar > sum(cases)) {
    stop("cstar should be at least 1 and less than or equal to the sum(cases)")
  }
}

#' Check longlat argument
#'
#' @param longlat A logical value indicating whether longlat distance should be used (TRUE)
#' @noRd
arg_check_longlat = function(longlat) {
  if (length(longlat) != 1) {
    stop("length(longlat) != 1")
  }
  if (!is.logical(longlat)) {
    stop("longlat should be a logical value")
  }
}

#' Check alpha argument
#'
#' @param alpha Signifcance level (single value > 0 and <= 1)
#' @noRd
arg_check_alpha = function(alpha) {
  if (length(alpha) != 1 || !is.numeric(alpha)) {
    stop("alpha should be a numeric vector of length 1")
  }
  if (alpha < 0 || alpha > 1) {
    stop("alpha should be a value [0, 1)")
  }
}

#' Check noc argument
#'
#' @param noc Logical value. Should only non-overlapping clusters be returned.
#' @noRd
arg_check_noc = function(noc) {
  if (length(noc) != 1) {
    stop("length(noc) != 1")
  }
  if (!is.logical(noc)) {
    stop("noc should be a logical value")
  }
}

#' Check ex argument
#'
#' @param ex A vector of expected counts
#' @param N nrow(coords)
#' @noRd
arg_check_ex = function(ex, N) {
  if (length(ex) != N) {
    stop("length(ex) != nrow(coords)")
  }
  if (!is.numeric(ex)) {
    stop("ex should be numeric values")
  }
  if (!is.vector(ex)) {
    stop("ex should be a vector")
  }
}

#' Check modified argument
#'
#' @param modified A logicical value. should bn.test be modified (TRUE)
#' @return NULL
#' @noRd
arg_check_modified = function(modified) {
  if (length(modified) != 1) {
    stop("length(modified) != 1")
  }
  if (!is.logical(modified)) {
    stop("modified should be a logical value")
  }
}

#' Check tobs argument
#'
#' For smerc_cluster function
#'
#' @param tobs Vector of observed test statistics
#' @return NULL
#' @noRd
arg_check_tobs = function(tobs) {
  if (!is.numeric(tobs)) {
    stop("tobs should be numeric values")
  }
  if (!is.vector(tobs)) {
    stop("tobs should be a vector")
  }
  if (min(tobs) < 0) {
    stop("tobs values must be >= 0")
  }
}

#' Check zones argument
#'
#' For smerc_cluster function
#'
#' @param zones A list of zones
#' @param N Number of tobs
#' @return NULL
#' @noRd
arg_check_zones = function(zones, N) {
  if (length(zones) != N) {
    stop("length(zones) != length(tobs)")
  }
  if (!list(zones)) {
    stop("zones should be a list")
  }
}

#' Check pvalue argument
#'
#' For smerc_cluster function
#'
#' @param pvalue Vector of p-values
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_pvalue = function(pvalue, N) {
  if (length(pvalue) != N) {
    stop("length(pvalue) != length(pvalue)")
  }
  if (!is.numeric(pvalue)) {
    stop("pvalue should be numeric values")
  }
  if (!is.vector(pvalue)) {
    stop("pvalue should be a vector")
  }
  if (min(pvalue) < 0 | max(pvalue) > 1) {
    stop("pvalue must have values in [0, 1]")
  }
}

arg_check_d = function(d, N) {
  if (!is.matrix(d)) {
    stop("d must be a matrix")
  }
  if (nrow(d) != N | ncol(d) != N) {
    stop("d must be square nrow(d) = nrow(coords)")
  }
}

#' Check method argument
#'
#' @param method A single character vector specifying
#' the method name.
#' @return NULL
#' @noRd
arg_check_method = function(method) {
  if (length(method) != 1) {
    stop("method must be a vector of length 1")
  }
  if (!is.vector(method)) {
    stop("method must be a vector of length 1")
  }
  if (!is.character(method)) {
    stop("method must be a character vector")
  }
}

#' Check method argument
#'
#' @param rel_param A list with relevant paramameters for a method
#' the method name.
#' @return NULL
#' @noRd
arg_check_rel_param = function(rel_param) {
  if (!is.list(rel_param)) {
    stop("rel_param must be a vector")
  }
}

#' Check shape_all argument
#'
#' For smerc_cluster function
#'
#' @param shape_all Vector of shapes
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_shape_all = function(shape_all, N) {
  if (length(shape_all) != N) {
    stop("length(shapes_all) != length(tobs)")
  }
  if (!is.numeric(shape_all)) {
    stop("shape_all should be numeric values")
  }
  if (!is.vector(shape_all)) {
    stop("shape_all should be a vector")
  }
  if (min(shape_all) < 1) {
    stop("All shapes must be >= 1")
  }
}

#' Check angle_all argument
#'
#' For smerc_cluster function
#'
#' @param angle_all Vector of shapes
#' @param N length(tobs)
#' @return NULL
#' @noRd
arg_check_angle_all = function(angle_all, N) {
  if (length(angle_all) != N) {
    stop("length(shapes_all) != length(tobs)")
  }
  if (!is.numeric(angle_all)) {
    stop("angle_all should be numeric values")
  }
  if (!is.vector(angle_all)) {
    stop("angle_all should be a vector")
  }
  if (min(angle_all) < 0 | max(angle_all) >= 180) {
    stop("All angles must be in (0, 180)")
  }
}

#' Check a argument
#'
#' @param a Penalty parameter for elliptic.test
#' @return NULL
#' @noRd
arg_check_a = function(a) {
  if (length(a) != 1) {
    stop("a must be a single value")
  }
  if (!is.numeric(a)) {
    stop("a must be a numeric value")
  }
  if (!is.vector(a)) {
    stop("a must be a vector")
  }
  if (a < 0) {
    stop("a must be >= 0")
  }
}

#' Check w argument
#'
#' @param w Spatial connectivity matrix
#' @param N nrow(coords)
#' @return NULL
#' @noRd
arg_check_w = function(w, N) {
  if (!is.matrix(w) | !is.data.frame(w)) {
    stop("w must be a matrix or data.frame")
  }
  if (nrow(w) != N | ncol(N)) {
    stop("w must be a square matrix with nrow(w) = nrow(coords)")
  }
}
