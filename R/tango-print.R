#' Print object of class \code{tango}.
#'
#' Print \code{tango} object
#'
#' @param x An object of class \code{tango}.
#' @param digits Number of significant digits to print.
#' @inheritDotParams base::print
#' @export
#' @examples
#' data(nydf)
#' coords = as.matrix(nydf[,c("x", "y")])
#' w = dweights(coords, kappa = 1)
#' results = tango.test(nydf$cases, nydf$pop, w, nsim = 49)
#' results
print.tango = function(x, ..., digits = 2) {
  message(paste(crayon::blue("method:"),
                crayon::magenta("Tango's index")))
  index = base::signif(x$tstat, digits = digits)
  gof = base::signif(x$gof, digits = digits)
  sa = base::signif(x$sa, digits = digits)
  tstat_chisq = base::signif(x$tstat.chisq, digits = digits)
  dfc = base::signif(x$dfc, digits = digits)
  pvalue_chisq = base::signif(x$pvalue.chisq, digits = digits)

  message(paste(crayon::blue("index:"),
                crayon::magenta(index)))
  message(paste(crayon::blue("goodness-of-fit component:"),
                crayon::magenta(gof)))
  message(paste(crayon::blue("spatial autocorrelation component:"),
                crayon::magenta(sa)))
  message(paste(crayon::blue("chi-square statistic:"),
                crayon::magenta(tstat_chisq)))
  message(paste(crayon::blue("chi-square df:"),
                crayon::magenta(dfc)))
  message(paste(crayon::blue("chi-square p-value:"),
                crayon::magenta(pvalue_chisq)))
  if (!is.null(x$pvalue.sim)) {
    pvalue_sim = base::signif(x$pvalue.sim, digits = digits)
    message(paste(crayon::blue("Monte Carlo p-value:"),
                  crayon::magenta(pvalue_sim)))
  }
}
