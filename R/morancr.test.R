#' Constant-risk Moran's I-based test
#'
#' \code{morancr.test} performs a test of clustering using the constant-risk
#' version of the Moran's I statistic proposed by Walter (1992) under the
#' constant risk hypothesis.
#'
#' @inheritParams flex.test
#' @param alternative a character string specifying the alternative hypothesis,
#'   must be one of "greater" (default), "two.sided", or "less". You can specify
#'   just the initial letter.
#' @return Returns a \code{smerc_similarity_test}.
#' @author Joshua French
#' @seealso \code{\link{morancr.stat}}
#' @export
#' @references Walter, S. D. (1992). The analysis of regional patterns in health
#'   data: I. Distributional considerations. American Journal of Epidemiology,
#'   136(6), 730-741.
#' @examples
#' data(nydf)
#' data(nyw)
#' morancr.test(cases = nydf$cases, pop = nydf$pop, w = nyw, nsim = 9)
morancr.test = function(cases, pop, w, ex = sum(cases) / sum(pop) * pop,
                         nsim = 499, alternative = "greater") {
  arg_check_morancr_stat(cases = cases, pop = pop, w = w, ex = ex)
  arg_check_nsim(nsim)
  alternative = match.arg(alternative, c("greater", "less", "two.sided"))
  tstat = morancr.stat(cases = cases, w = w, ex = ex)
  tsim = morancr.sim(nsim = nsim, cases = cases, w = w, ex = ex)
  pvalue = NA
  if (alternative == "two.sided") {
    pvalue = (sum(abs(tsim) >= abs(tstat)) + 1) / (nsim + 1)
  } else if (alternative == "greater") {
    pvalue = (sum(tsim >= tstat) + 1) / (nsim + 1)
  } else if (alternative == "less") {
    pvalue = (sum(tsim <= tstat) + 1) / (nsim + 1)
  }
  return(structure(list(test_statistic = tstat,
                        tsim = tsim,
                        pvalue = pvalue,
                        alternative = alternative,
                        nsim = nsim,
                        simdist = "multinomial"), class = "smerc_similarity_test"))
}
