\name{dsT_covariates}
\docType{data}
\alias{dsT_covariates}
\title{Environmental covariate suite of rasters for testing DSMART}
\description{A \code{RasterStack} of 20 environmental covariate rasters. These data are described in Odgers et al. (2014). The assembled covariates are surrogates for soil forming factors and are commonly referred to as scorpan factors, where scorpan is a mnemonic of the soil forming factors (s) soil, (c) climate, (o) organism, (r) relief, (p) parent material, (a) age, and (n) spatial position (McBratney et al. 2003). The \code{RasterStack} is a small subset from the Dalrymple Shire in Central Queensland, Australia. The subset area is a 7km by 7km square, and the rasters have a 30m pixel resolution. The coordinate reference system of these maps is GDA 1994 Australia Albers.}

\usage{data(dsT_covariates)}
\format{
\code{dsT_covariates} is a 20 layer \code{RasterStack}.}
\details{This is an example data set in order to familiarize users to the working and inputs required for operationalising DSMART which is implemented in the companion functions \code{dsmart} and \code{dsmartR}. }

\references{
\itemize{
\item McBratney, A.B., Mendonca Santos, M.L., Minasny, B., (2003) \href{http://dx.doi.org/10.1016/S0016-7061(03)00223-4}{On digital soil mapping}. Geoderma, 117: 3-52.
\item Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014) \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map units through resampled classification trees}. Geoderma, 214-215: 91-100.
}
}

\examples{
library(dsmart)
library(raster)

data(dsT_covariates)
names(dsT_covariates)
}
\keyword{datasets}
