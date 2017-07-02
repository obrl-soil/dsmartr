\name{dsT_composition}
\docType{data}
\alias{dsT_composition}
\title{Soil mapping unit compositions}
\description{A \code{data.frame} that tabulates the soil class compositions for each mapping unit from the \code{dsT_polyons} soil map. The nature of this data is a required input for the \code{dsmart} function. This data is for testing the companion \code{dsmart} and \code{dsmartR} functions.}

\usage{data(dsT_composition)}
\format{
\code{dsT_composition} is of class \code{data.frame}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}. It is composed of 4 columns. The first is the numeric identifier of a map unit, the second is the map unit character name, the third are the soil class codes and the last are the numeric proportions of each class to each soil mapping unit. All map units and their compositions are row-binded together.}

\references{
\itemize{
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(dsmart)

data(dsT_composition)
levels(as.factor(dsT_composition$poly))
levels(as.factor(dsT_composition$soil_class))

}
\keyword{datasets}
