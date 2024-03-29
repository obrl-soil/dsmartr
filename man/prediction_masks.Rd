% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dsmartr_helpers.R
\name{prediction_masks}
\alias{prediction_masks}
\title{generate masks from samples}
\usage{
prediction_masks(samples = NULL, covariates = NULL, tolerance = 0L, cpus = 1)
}
\arguments{
\item{samples}{List; POINT_sfc objects created by
\code{\link[dsmartr:iterate]{dsmartr::iterate()}} when option write_samples
= TRUE}

\item{covariates}{RasterStack or RasterBrick; environmental covariate
datasets.}

\item{tolerance}{Integer; the number of out-of-range covariates allowed.
Defaults to 0, that is, if any covariate is out of range on a given pixel,
that pixel will be masked.}

\item{cpus}{Integer; number of processors to use in parallel.}
}
\value{
a list of rasters
}
\description{
Generates a rasterlayer for each dsmart iteration highlighting areas where
one or more covariates has a value that is out of range of the sample used to
model soil distribution.
}
\details{
For each model iteration this function generates a raster that can
  be used to mask dsmartr predictions where covariate values are out of range
  of the sample of cells used in that iteration. The masks can also be
  combined to produce an overall mask, applicable to final products. This is
  experimental; the idea is to avoid making predictions in locations where
  the model won't perform well.
}
\examples{
\dontrun{
# run dsmartr::iterate() with the example data, then:
sample_list <- lapply(list.files(file.path(getwd(), 'iterations', 'samples'),
pattern = '\\\\.gpkg', full.names = TRUE), read_sf)
masks <- prediction_masks(samples = sample_list,
                          covariates = heronvale_covariates,
                          cpus = max(1, (parallel::detectCores() - 1)))

# less strict: up to 3 out of range covariates are allowed
masks <- prediction_masks(samples = sample_list,
                          covariates = heronvale_covariates,
                          tolerance = 3,
                          cpus = max(1, (parallel::detectCores() - 1)))

# a 'final product' mask, where outputs are masked only if the mode of that
# cell is 'do not predict'.
all_masks <- stack(masks)
modal_mask <- calc(all_masks, function(cell) ifelse(raster::modal(cell) == 1,
NA, 0))
# the above could be applied to the most-likely soils map after running
# dsmartr::most_likely:
masked_m1 <- most_likely_soil[[1]][[1]] + modal_mask
}
}
