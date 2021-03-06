% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsmartr_unstack.R
\name{most_likely}
\alias{most_likely}
\title{Extract most-likely-soil maps}
\usage{
most_likely(dsmart_preds = NULL, dsmart_probs = NULL, n_maps = 2L)
}
\arguments{
\item{dsmart_preds}{RasterBrick; output \code{dsmartr_predictions} from
\code{\link[dsmartr:collate]{dsmartr::collate()}}.}

\item{dsmart_probs}{RasterBrick; output \code{dsmartr_probabilities} from
\code{\link[dsmartr:collate]{dsmartr::collate()}}; optional.}

\item{n_maps}{integer; the number of most-probable maps desired.}
}
\value{
A list containing \code{n_maps} most-likely-soils maps, from most to
  least probable. Optionally; \code{n_maps} probability surfaces associated
  with \code{most_likely_maps} follow.
All outputs are written to disk as GeoTIFFs.
}
\description{
This function takes the outputs of \code{\link{collate}} and extracts the top
\code{n} most-likely soil maps. Optionally, their respective probability
surfaces are also extracted.
}
\examples{
\dontrun{
# run collate() example code, then:
most_likely_soil <-
 most_likely(dsmart_preds = collated[['dsmartr_predictions']],
             dsmart_probs = collated[['dsmartr_probabilities']],
             n_maps = 2)
}
}
