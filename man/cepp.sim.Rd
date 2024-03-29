% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cepp.sim.R
\name{cepp.sim}
\alias{cepp.sim}
\title{Perform \code{cepp.test} on simulated data}
\usage{
cepp.sim(nsim = 1, nn, ty, ex, wts, simdist = "multinomial")
}
\arguments{
\item{nsim}{A positive integer indicating the number of
simulations to perform.}

\item{nn}{A list of nearest neighbors produced by \code{\link{casewin}}.}

\item{ty}{The total number of cases in the study area.}

\item{ex}{The expected number of cases for each region.
The default is calculated under the constant risk
hypothesis.}

\item{wts}{A list that has the weights associated with each
region of each element of \code{nn}.}

\item{simdist}{A character string indicating whether the
simulated data should come from a \code{"multinomial"}
or \code{"poisson"} distribution.  The default is
\code{"multinomial"}, which fixes the total number of
cases observed in each simulated data set.}
}
\value{
A vector with the maximum test statistic for each
  simulated data set.
}
\description{
\code{cepp.sim} efficiently performs
\code{\link{cepp.test}} on a simulated data set.  The
function is meant to be used internally by the
\code{\link{cepp.test}} function, but is informative for
better understanding the implementation of the test.
}
\examples{
data(nydf)
coords <- with(nydf, cbind(longitude, latitude))
d <- gedist(as.matrix(coords), longlat = TRUE)
nn <- casewin(d, cases = nydf$pop, cstar = 15000)
cases <- floor(nydf$cases)
ty <- sum(cases)
ex <- ty / sum(nydf$pop) * nydf$pop
# find smallest windows with at least n* pop
nstar <- 1000
nn <- casewin(d, cases = nydf$pop, cstar = nstar)
# determine ts
wts <- cepp.weights(nn, nydf$pop, nstar)
tsim <- cepp.sim(1, nn = nn, ty = ty, ex = ex, wts = wts)
}
