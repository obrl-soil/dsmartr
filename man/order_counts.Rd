% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsmartr_collate.R
\name{order_counts}
\alias{order_counts}
\title{Order soil class counts}
\usage{
order_counts(input = NULL, n_classes = NULL)
}
\arguments{
\item{input}{A vector of integers n model runs long, where the integers refer
to counts of soil classes.}

\item{n_classes}{Integer; the total number of soil classes available to
predict.}
}
\value{
An atomic vector of probabilities, or a vector of NA values of the
  same length. Ties are shuffled at random, e.g. for 20 model runs, e.g. 0 5
  8 0 becomes either 3 2 1 4 or 3 2 4 1. These values should correspond to a
  lookup table ID column, so they can be linked to a soil class name.
}
\description{
Rearranges soil class predictions into descending order, handling ties
sensibly.
}
\note{
This is a helper function and not widely applicable. It is expected to
  be used in concert with \code{\link[raster:calc]{raster::calc()}}, applied
  to a stack of tabulated model realisations produced by
  \code{\link[dsmartr:collate]{dsmartr::collate()}}.
}
\examples{
\dontrun{
# 100 model runs, from a cell in the middle of the heronvale demo dataset
cell_19291 <- c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9,
15, 34, 9, 24, 34, 24, 34, 21, 24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15,
21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24, 12, 21, 15, 15, 29, 12, 7, 21,
12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29, 7, 3, 34, 9, 9,
33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
15, 7)
preds_19291 <- count_predictions(cell_19291)
ordered_preds <- order_counts(preds_19291)
}

}
\keyword{internal}
