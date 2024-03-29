% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nclusters.R
\name{nclusters}
\alias{nclusters}
\title{Number of clusters}
\usage{
nclusters(x)
}
\arguments{
\item{x}{A \code{smerc_cluster object}}
}
\value{
A non-negative integer.
}
\description{
\code{nclusters} returns the number of clusters
identified in a \code{smerc_cluster} object.
}
\examples{
data(nydf)
coords <- with(nydf, cbind(longitude, latitude))
out <- scan.test(
  coords = coords, cases = floor(nydf$cases),
  pop = nydf$pop, nsim = 19,
  alpha = 0.3, longlat = TRUE
)
nclusters(out)
}
