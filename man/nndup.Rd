% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nndup.R
\name{nndup}
\alias{nndup}
\title{Determine duplicates in nearest neighbor list}
\usage{
nndup(nn, N = max(unlist(nn)))
}
\arguments{
\item{nn}{A list of nearest neighbors.}

\item{N}{The largest value in \code{nn}.}
}
\value{
A logical vector of indicating duplicate indices.
}
\description{
\code{nndup} determines the indices of duplicated
elements for a nearest neighbors list created by a
function such as \code{\link{nnpop}} or
\code{\link{knn}}.  The indices are related to the list
returned by \code{\link{nn2zones}}.
}
\examples{
nn <- list(1:3, c(2:1, 4))
nndup(nn, 4)
}
