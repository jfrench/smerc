#' Print object of class \code{smerc_cluster}.
#'
#' Print \code{scan} object
#'
#' @param x An object of class \code{smerc_cluster}.
#' @inheritDotParams base::print
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' out
print.smerc_cluster = function(x, ...) {
  message(paste(crayon::blue("method:"),
                crayon::magenta(x$method)))
  rel_param = x$rel_param
  if (!is.null(rel_param$type)) {
    message(paste(crayon::blue("statistic:"),
                  crayon::magenta(rel_param$type)))
  }
  if (!is.null(rel_param$simdist)) {
    message(paste(crayon::blue("simulation:"),
                  crayon::magenta(rel_param$simdist)))
  }
  if (!is.null(rel_param$nsim)) {
    message(paste(crayon::blue("realizations:"),
                  crayon::magenta(rel_param$nsim)))
  }
  if (!is.null(rel_param$longlat)) {
    message(paste(crayon::blue("distance:"),
            ifelse(rel_param$longlat,
                   crayon::magenta("great circle"),
                   crayon::magenta("euclidean"))
            ))
  }
  if (!is.null(rel_param$ubpop)) {
    message(paste(crayon::blue("population upperbound: "),
                  crayon::magenta(rel_param$ubpop * 100, "%",
                                  sep = "")))
  }
  if (!is.null(rel_param$min.cases)) {
    message(paste(crayon::blue("minimum cases:"),
                  crayon::magenta(rel_param$min.cases)))
  }
  if (!is.null(rel_param$alpha)) {
    message(paste(crayon::blue("significance level:"),
                  crayon::magenta(rel_param$alpha)))
  }
}
