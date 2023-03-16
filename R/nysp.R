#' @name nysp
#' @title \code{SpatialPolygonsDataFrame} for New York leukemia data.
#' @description The number of incident leukemia cases from 1978-1982 per census tract for an 8-county region
#' of upstate New York.
#'
#' This is the same data as in \code{nydf} in a different format.
#'
#' Note that the coordinates in the polygons have been projected to a different
#' coordinate system (UTM, zone 18) compared to \code{nydf}, but the order of the
#' regions/polygons is the same as in \code{nydf}.
#'
#' @format A \code{SpatialPolygonsDataFrame} with 281 rows and 17 columns:
#' \describe{
#'  \item{areaname}{The name of the region.}
#'  \item{areakey}{Census tract id.}
#'  \item{xorig}{x-coordinate associated with the centroid of each region on the ORIGINAL scale.}
#'  \item{yirug}{y-coordinate associated with the centroid of each region on the ORIGINAL scale.}
#'  \item{pop8}{The population (1980 census) of the region.}
#'  \item{tractcas}{The number of leukemia cases between 1978-1982, rounded to two decimal places.}
#'  \item{propcas}{The proportion of cases relative to population.}
#'  \item{pctownhome}{The percentage of homeowners in the tract.}
#'  \item{pctage65p}{The percentage of residents aged 65 or older.}
#'  \item{Z}{A transformation of exposure to TCE, specifically log(1000(TCE + 1)/pop8).}
#'  \item{avgidist}{Average inverse distance to the nearest TCE site.}
#'  \item{pexpossure}{}
#'  \item{cases}{The number of leukemia cases between 1978-1982.}
#'  \item{xm}{A shifted version of \code{x}.}
#'  \item{ym}{A shifted version of \code{y}.}
#'  \item{x}{x-coordinate associated with the centroid of each region.}
#'  \item{y}{y-coordinate associated with the centroid of each region.}
#' }
#' @source Waller, L.A. and Gotway, C.A. (2005).
#' Applied Spatial Statistics for Public Health Data.
#' Hoboken, NJ: Wiley.
#'
#' @docType data
#' @usage data(nysp)
#'
#' @source Bivand, R. S., Pebesma, E. J., Gomez-Rubio, V.,
#' and Pebesma, E. J. (2013). Applied Spatial Data Analysis
#' with R, 2nd edition. New York: Springer.
NULL
