% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smerc_cluster-plot.R
\name{plot.smerc_cluster}
\alias{plot.smerc_cluster}
\title{Plot object of class \code{smerc_cluster}.}
\usage{
\method{plot}{smerc_cluster}(
  x,
  ...,
  idx = seq_along(x$clusters),
  nclusters = NULL,
  ccol = NULL,
  cpch = NULL,
  add = FALSE,
  usemap = FALSE,
  mapargs = list()
)
}
\arguments{
\item{x}{An object of class scan to be plotted.}

\item{...}{Additional graphical parameters passed to the
\code{plot} function.}

\item{idx}{An index vector indicating the elements of
\code{x$clusters} to print information for. The default
is all clusters.}

\item{nclusters}{Number of clusters to plot. Deprecated. Use \code{idx}.}

\item{ccol}{Fill color of the plotted points.  Default is
\code{grDevices::hcl.colors(nclusters, palette = "viridis")}.}

\item{cpch}{Plotting character to use for points in each
cluster.  Default is NULL, indicating pch = 20 for the
most likely cluster and then pch = 2, 3, .., up to the
remaining number of clusters.}

\item{add}{A logical indicating whether results should be
drawn on existing map.}

\item{usemap}{Logical indicating whether the maps::map
function should be used to create a plot background for
the coordinates.  Default is \code{FALSE}.  Use
\code{TRUE} if you have longitude/latitude coordinates.}

\item{mapargs}{A list of arguments for the map function.}
}
\description{
Plot clusters (the centroids of the regions in each
cluster) in different colors.  The most likely cluster is
plotted with solid red circles by default.  Points not in
a cluster are black open circles.  The other cluster
points are plotted with different symbols and colors.
}
\examples{
data(nydf)
coords <- with(nydf, cbind(longitude, latitude))
out <- scan.test(
  coords = coords, cases = floor(nydf$cases),
  pop = nydf$pop, nsim = 0,
  longlat = TRUE, alpha = 1
)
# plot only 3 most likely clusters
plot(out, idx = 1:3)
## plot output for new york state
# specify desired argument values
mapargs <- list(
  database = "county", region = "new york",
  xlim = range(out$coords[, 1]),
  ylim = range(out$coords[, 2])
)
# needed for "county" database (unless you execute library(maps))
data(countyMapEnv, package = "maps")
# plot only the 1st and 3rd clusters
plot(out, idx = 1:3, usemap = TRUE, mapargs = mapargs)
}
