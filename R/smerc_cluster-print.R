#' Print object of class \code{smerc_cluster}.
#'
#' Print \code{scan} object
#'
#' @param x An object of class \code{smerc_cluster}.
#' @inheritDotParams base::print
#' @param extra A logical value. Default is \code{FALSE}.
#' \code{TRUE} indicates that extra information should be
#' printed.
#' @export
#' @examples
#' data(nydf)
#' coords = with(nydf, cbind(x, y))
#' out = scan.test(coords = coords, cases = floor(nydf$cases),
#'                 pop = nydf$pop, nsim = 49,
#'                 longlat = TRUE, alpha = 0.12)
#' out
print.smerc_cluster = function(x, ..., extra = FALSE) {
  if (length(extra) != 1 | !is.logical(extra)) {
    stop("extra must be a single logical")
  }
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
  if (!is.null(rel_param$ubpop)) {
    message(paste(crayon::blue("population upperbound: "),
                  crayon::magenta(rel_param$ubpop * 100, "%",
                                  sep = "")))
  }
  if (!is.null(rel_param$min.cases)) {
    message(paste(crayon::blue("minimum cases:"),
                  crayon::magenta(rel_param$min.cases)))
  }
  if (!is.null(rel_param$cstar)) {
    message(paste(crayon::blue("case radius:"),
                  crayon::magenta(rel_param$cstar)))
  }
  if (!is.null(rel_param$modified)) {
    message(paste(crayon::blue("modified p-value:"),
                  crayon::magenta(rel_param$modified)))
  }
  if (extra) {
    message(paste(crayon::blue("significance level:"),
                  crayon::magenta(x$alpha)))
    dtype = ifelse(x$longlat, "great circle", "euclidean")
    message(paste(crayon::blue("distance:"),
                  crayon::magenta(dtype)))
  }
}
