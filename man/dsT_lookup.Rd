\name{dsT_lookup}
\docType{data}
\alias{dsT_lookup}
\title{Lookup table for soil classes (DSMART test).}
\description{A \code{data.frame} that is automatically generated from the \code{dsmart} function that acts as a link between the character-based soil class labeling and its numerical counterpart i.e. a lookup table.  Using the other data objects in this package, namely \code{dsT_covariates}, \code{dsT_composition}, and \code{dsT_polygons}, the \code{dsmart} algorithm was set to run and generate outputs for 5 resampling and C5 modeling realizations. The lookup table was produced as a result of these parametisations. }

\usage{data(dsT_lookup)}
\format{
\code{dsT_lookup} is of class \code{data.frame}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}. It is composed of 2 columns. The first is the character identifier of a soil class, the second is its numerical counterpart of which is acutally stored for each soil class raster.}

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.}}

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
