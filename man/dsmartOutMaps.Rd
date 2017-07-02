\name{dsmartOutMaps}
\docType{data}
\alias{dsmartOutMaps}
\title{Companion data set for \code{dsmartR}: Outputs from \code{dsmart} function.}
\description{A \code{RasterStack} of 5 realizations from the \code{dsmart} function. Using the other data objects in this package, namely \code{dsT_covariates}, \code{dsT_composition}, and \code{dsT_polygons}, the \code{dsmart} algorithm was set to run and generate outputs for 5 resampling and C5 modeling realizations. The input data to produce these rasters are described in Odgers et al. (2014). The data are a small subset of the totally area described, which is the Dalrymple Shire in Central Queensland, Australia. The subset area is a 7km by 7km square, and the rasters have a 30m pixel resolution. The coordinate reference system of these maps is GDA 1994 Australia Albers.}

\usage{data(dsmartOutMaps)}
\format{
\code{dsmartOutMaps} is a 5 layer \code{RasterStack} of raster outputs from the \code{dsmart} algorithm.
}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}.}

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(dsmart)
library(raster)
library(rasterVis)

data(dsmartOutMaps)
data(dsT_lookup)

#Plot first realisation
map1<- dsmartOutMaps[[1]]
levelplot(map1)

#Not run
#To use the character based soil class labels for mapping legend
#map1<- as.factor(map1)
# Add soil class labels
#rat<- levels(map1)[[1]]
#rat[["soil class"]]<- dsT_lookup[which(dsT_lookup[,2] %in% c(as.matrix(rat))),1]
#levels(map1)<- rat
#levelplot(map1)

}
\keyword{datasets}
