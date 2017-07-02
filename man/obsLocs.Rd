\name{obsLocs}
\docType{data}
\alias{obsLocs}
\title{A collection of soil class observations}
\description{A \code{data.frame} containing point observations of soil classes. This data is for demonstration of the usage of incorporating observations within the soil map disaggregation procedure.}

\usage{data(obsLocs)}
\format{
\code{obsLocs} is of class \code{data.frame}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}. It is composed of 3 columns. Columns 1 and 2 are spatial coordinates. Column 3 are soil class labels.}

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(dsmart)
library(raster)

#load the relevent data

#Covariates
data(dsT_covariates)

#Polygons
data(dsT_polygons)

#Map unit compositions
data(dsT_composition)

##Not Run
#Run dsmart (with 15 samples per polygon, 5 C5 model realisations, using 1 compute node)
# Incorporate addtional point observations
#data(obsLocs)
#dsmart(covariates = dsT_covariates, polygons = dsT_polygons, composition = dsT_composition, obsdat=obsLocs, n=15, reals = 5, cpus=1)

}
\keyword{datasets}
