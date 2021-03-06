% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsmartr_collate.R
\name{calc_probabilities}
\alias{calc_probabilities}
\title{Soil class probabilities}
\usage{
calc_probabilities(input = NULL, n_iterations = NULL)
}
\arguments{
\item{input}{A vector of integers n soil classes long, where the integers
  refer to counts of soil
classes.}

\item{n_iterations}{Integer; the total number of model runs.}
}
\value{
An atomic vector of probabilities, or a vector of NA values of the
  same length. e.g. for 20 model runs, 0 5 8 0 becomes 0 0.25 0.4 0.
}
\description{
Converts the number of times a soil class was predicted to a probability by
normalising against the total number of model runs.
}
\note{
This is a helper function and not widely applicable. It is expected to
  be used in concert with \code{\link[raster:calc]{calc}}, applied to a stack
  of tabulated class counts produced by
  \code{\link[dsmartr:collate]{collate}}.
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
probs_19291 <- calc_probabilities(preds_19291)
}

}
\keyword{internal}
