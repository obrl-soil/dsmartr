\name{dsT_polygons}
\docType{data}
\alias{dsT_polygons}
\title{Soil polygon map}
\description{A \code{SpatialPolygonsDataFrame} of a soil map showing soil mapping units. This map is essentially a common shapefile that has been read into the \code{R} environment. This data is for testing the companion \code{dsmart} and \code{dsmartR} functions and is a subset of the soil map described and disaggregated by Odgers et al. (2014). The polygons represent soil mapping units. The attribute table contains information regarding the soil class compositions and their proportions within each soil mapping unit.}

\usage{data(dsT_polygons)}
\format{
\code{dsT_polygons} is of class \code{SpatialPolygonsDataFrame}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}. }

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(dsmart)
library(sp)

data(dsT_polygons)
#make a plot
plot(dsT_polygons)
}
\keyword{datasets}
