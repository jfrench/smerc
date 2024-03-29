% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rflex.sim.R
\name{rflex.sim}
\alias{rflex.sim}
\title{Perform \code{rflex.test} on simualated data}
\usage{
rflex.sim(
  nsim = 1,
  nn,
  w,
  ex,
  alpha1 = 0.2,
  type = "poisson",
  pop = NULL,
  cl = NULL
)
}
\arguments{
\item{nsim}{A positive integer indicating the number of
simulations to perform.}

\item{nn}{A matrix of the k nearest neighbors for the
regions described by \code{w}.}

\item{w}{A binary spatial adjacency matrix for the
regions.}

\item{ex}{The expected number of cases for each region.
The default is calculated under the constant risk
hypothesis.}

\item{alpha1}{The middle p-value threshold.}

\item{type}{The type of scan statistic to compute. The
default is \code{"poisson"}. The other choice
is \code{"binomial"}.}

\item{pop}{The population size associated with each
region.}

\item{cl}{
A cluster object created by \code{\link{makeCluster}},
or an integer to indicate number of child-processes
(integer values are ignored on Windows) for parallel evaluations
(see Details on performance).
It can also be \code{"future"} to use a future backend (see Details),
\code{NULL} (default) refers to sequential evaluation.
}
}
\value{
A vector with the maximum test statistic for each
  simulated data set.
}
\description{
\code{rflex.sim} efficiently performs
\code{\link{rflex.test}} on a simulated data set.  The
function is meant to be used internally by the
\code{\link{rflex.test}} function, but is informative for
better understanding the implementation of the test.
}
\examples{
data(nydf)
data(nyw)
# determine knn
coords <- with(nydf, cbind(longitude, latitude))
nn <- knn(coords, longlat = TRUE, k = 50)
# determine expected number of cases in each region
cases <- floor(nydf$cases)
pop <- nydf$pop
ex <- pop * sum(cases) / sum(pop)
tsim <- rflex.sim(nsim = 5, nn = nn, w = nyw, ex = ex)
}
