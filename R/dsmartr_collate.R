#' Soil class prediction counts
#'
#' Counts the number of times a given soil class was predicted on a pixel
#' @keywords internal
#' @param input A vector of integers n model runs long, where the integers refer
#'   to soil classes.
#' @param n_classes Integer; the total number of soil classes available to
#'   predict.
#' @return An atomic vector of prediction counts, or a vector of NA values of
#'   the same length. e.g. 0 5 8 0 == soil 1 was not predicted, soil 2 was
#'   predicted 5 times, soil 3 was predicted 8 times, soil 4 was not predicted.
#' @note This is a helper function and not widely applicable. It is expected to
#'   be used in concert with \code{\link[raster:calc]{calc}}, applied to a stack
#'   of model realisations produced by \code{\link[dsmartr:iterate]{iterate}}.
#' @examples \dontrun{
#' # 100 model runs, from a cell in the middle of the heronvale demo dataset
#' cell_19291 <- c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9,
#' 15, 34, 9, 24, 34, 24, 34, 21, 24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15,
#' 21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24, 12, 21, 15, 15, 29, 12, 7, 21,
#' 12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29, 7, 3, 34, 9, 9,
#' 33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
#' 15, 7)
#' count_predictions(cell_19291)
#' }
#'
count_predictions <- function(input = NULL, n_classes = NULL) {
  # NB is.na(sum((input))) is faster but does not allow for cases where some
  # runs don't make a prediction on a pixel and others do. This can happen with
  # covariates of varying extent (e.g. along a coastline), but only in cases
  # where some special options in C50::C50Control have been used, e.g. winnow =
  # TRUE). Bear in mind that in these cases, probabilities are still normalised
  # against total model runs.

  if (all(is.na(input))) {
    rep(NA_integer_, n_classes)
  } else  {
    tabulate(input, nbins = n_classes)
  }
}

#' Soil class probabilities
#'
#' Converts the number of times a soil class was predicted to a probability by
#' normalising against the total number of model runs.
#' @keywords internal
#' @param input A vector of integers n soil classes long, where the integers
#'   refer to counts of soil
#' classes.
#' @param n_iterations Integer; the total number of model runs.
#' @return An atomic vector of probabilities, or a vector of NA values of the
#'   same length. e.g. for 20 model runs, 0 5 8 0 becomes 0 0.25 0.4 0.
#' @note This is a helper function and not widely applicable. It is expected to
#'   be used in concert with \code{\link[raster:calc]{calc}}, applied to a stack
#'   of tabulated class counts produced by
#'   \code{\link[dsmartr:collate]{collate}}.
#' @examples \dontrun{
#' # 100 model runs, from a cell in the middle of the heronvale demo dataset
#' cell_19291 <- c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9,
#' 15, 34, 9, 24, 34, 24, 34, 21, 24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15,
#' 21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24, 12, 21, 15, 15, 29, 12, 7, 21,
#' 12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29, 7, 3, 34, 9, 9,
#' 33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
#' 15, 7)
#' preds_19291 <- count_predictions(cell_19291)
#' probs_19291 <- calc_probabilities(preds_19291)
#' }
#'
calc_probabilities <- function(input = NULL, n_iterations = NULL) {
  round((input / n_iterations), 3)
  }

#' Order soil class counts
#'
#' Rearranges soil class predictions into descending order, handling ties
#' sensibly.
#' @keywords internal
#' @param input A vector of integers n model runs long, where the integers refer
#'   to counts of soil classes.
#' @param n_classes Integer; the total number of soil classes available to
#'   predict.
#' @return An atomic vector of probabilities, or a vector of NA values of the
#'   same length. Ties are shuffled at random, e.g. for 20 model runs, e.g. 0 5
#'   8 0 becomes either 3 2 1 4 or 3 2 4 1. These values should correspond to a
#'   lookup table ID column, so they can be linked to a soil class name.
#' @note This is a helper function and not widely applicable. It is expected to
#'   be used in concert with \code{\link[raster:calc]{raster::calc()}}, applied
#'   to a stack of tabulated model realisations produced by
#'   \code{\link[dsmartr:collate]{dsmartr::collate()}}.
#' @examples \dontrun{
#' # 100 model runs, from a cell in the middle of the heronvale demo dataset
#' cell_19291 <- c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9,
#' 15, 34, 9, 24, 34, 24, 34, 21, 24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15,
#' 21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24, 12, 21, 15, 15, 29, 12, 7, 21,
#' 12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29, 7, 3, 34, 9, 9,
#' 33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
#' 15, 7)
#' preds_19291 <- count_predictions(cell_19291)
#' ordered_preds <- order_counts(preds_19291)
#' }
#'
order_counts <- function(input = NULL, n_classes = NULL) {
  if (is.na(sum(input))) {
    rep(NA_integer_, n_classes)
    } else {
      order(input, sample(input, size = n_classes, replace = FALSE),
            decreasing = TRUE, na.last = TRUE)
    }
}

#' Sort soil class probabilities
#'
#' Rearranges soil class probabilities into descending order.
#' @keywords internal
#' @param input A vector of integers n model runs long, where the integers refer
#'   to soil class occurrence probability.
#' @param n_classes Integer; the total number of soil classes available to
#'   predict.
#' @return An atomic vector of probabilities, or a vector of NA values of the
#'   same length. Ties are shuffled at random, e.g. for 20 model runs, e.g. 0
#'   0.25 0.4 0 becomes 0.4 0.25 0 0.
#' @note This is a helper function and not widely applicable. It is expected to
#'   be used in concert with \code{\link[raster:calc]{raster::calc()}}, applied
#'   to a stack of tabulated model realisations produced by
#'   \code{\link[dsmartr:collate]{dsamrtr::collate()}}.
#' @examples \dontrun{
#' # 100 model runs, from a cell in the middle of the heronvale demo dataset
#' cell_19291 <- c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9,
#' 15, 34, 9, 24, 34, 24, 34, 21, 24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15,
#' 21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24, 12, 21, 15, 15, 29, 12, 7, 21,
#' 12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29, 7, 3, 34, 9, 9,
#' 33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
#' 15, 7)
#' preds_19291 <- count_predictions(cell_19291)
#' probs_19291 <- calc_probabilities(preds_19291)
#' ordered_preds <- order_counts(preds_19291)
#' ordered_probs <- sort_probabilities(probs_19291)
#' }
#'
sort_probabilities <- function(input = NULL, n_classes = NULL) {
  if (is.na(sum(input))) {
    # nb supplying n_classes is faster than length(input)
    rep(NA_real_, n_classes)
    } else {
      sort(input, decreasing = TRUE, na.last = TRUE)
    }
}

#' Collate dsmartr iterations
#'
#' Processes the outputs of \code{\link[dsmartr:iterate]{iterate}}
#' @param iteration_stack RasterStack or Brick; output \code{iteration_maps}
#'   from \code{\link{iterate}}.
#' @param lookup Data Frame; contains raster values and corresponding soil class
#'   labels. Example:
#'   \preformatted{'data.frame':	38 obs. of  2 variables:
#'   $ VALUE: int  1 2 3 4 5 6 7 8 9 10 ...
#'   $ CLASS: chr  "Ad" "An" "Bb" "Bh" ...}
#' @param cpus Integer; number of processors to use in parallel.
#' @return A list of four RasterStacks:
#' \itemize{
#'   \item{\code{dsmartr_predictions}: RasterStack, soil class prediction maps
#'   in order of most to least-probable}
#'   \item{\code{dsmartr_probabilties}: RasterStack, a probability surface for
#'   each layer of \code{dsmartr_predictions}.}
#'   \item{\code{tallied_predictions}: RasterStack containing per-pixel tallies
#'   of soil class occurrence.}
#'   \item{\code{tallied_probabilities}: RasterStack containing per-pixel
#'   probabilities of soil class occurrence}.}
#'
#' All outputs are written to disk as multiband GeoTIFFs before being assigned
#' to the global environment.
#' @note This function can generate very large R temporary files. A rough
#'   minimum disk space requirement is \code{n} iterations * \code{n} covariate
#'   cells * 8 bytes * 2.
#' @examples \dontrun{
#' # run iterate() with the example data, then:
#' LUT <- levels(iteration_maps[[1]])[[1]]
#' collated <- collate(iteration_stack = iteration_maps, lookup = LUT,
#'  cpus = max(1, (parallel::detectCores() - 1)))
#'  }
#' @importFrom raster beginCluster calc clusterR endCluster getCluster nlayers
#'   stack
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
collate <- function(iteration_stack = NULL,
                    lookup          = NULL,
                    cpus            = 1) {
  message(paste0(Sys.time(), ': dsmartr collation in progress...'))

  if (!dir.exists(file.path(getwd(), 'tallies'))) {
    dir.create(file.path(getwd(), 'tallies'), showWarnings = FALSE)
  }
  strs         <- file.path(getwd(), 'tallies')
  soil_classes <- dim(lookup)[1]
  n_iterations <- nlayers(iteration_stack)

  beginCluster(n = cpus)
  assign('soil_classes', soil_classes, envir = parent.frame())
  assign('n_iterations', n_iterations, envir = parent.frame())
  # 1. tally up how many times each soil class was predicted
  counts <- clusterR(x         = iteration_stack,
                     fun       = calc,
                     args      = list(fun = function(cell) {
                       count_predictions(input = cell, n_classes = soil_classes)
                       }),
                     filename  = file.path(strs, 'tallied_predictions.tif'),
                     datatype  = 'INT2S',
                     format    = 'GTiff',
                     NAflag    = -9999,
                     overwrite = TRUE)

  message(paste0(Sys.time(),
                 ': Predictions tallied. Calculating probabilities...'))

  # 2. express the tallies as probablities
  probs <- clusterR(x         = counts,
                    fun       = calc,
                    args      = list(fun = function(cell) {
                      calc_probabilities(input = cell,
                                         n_iterations = n_iterations)
                      }),
                    filename  = file.path(strs, 'tallied_probabilities.tif'),
                    datatype  = 'FLT4S',
                    format    = 'GTiff',
                    NAflag    = -9999,
                    overwrite = TRUE)

  message(paste0(Sys.time(),
                 ': Probabilities calculated. Ordering predictions...'))

  # 3. get the soil classes corresponding to counts
  ordered <- clusterR(x         = counts,
                      fun       = calc,
                      args      = list(fun = function(cell) {
                        order_counts(input = cell, n_classes = soil_classes)
                      }),
                      filename  = file.path(strs, 'dsmartr_predictions.tif'),
                      datatype  = 'INT2S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)

  message(paste0(Sys.time(),
                 ': Predictions ordered. Ordering probabilities...'))

  # 4. sort probs by most to least probable
  sorted  <- clusterR(x         = probs,
                      fun       = calc,
                      args      = list(fun = function(cell) {
                        sort_probabilities(input = cell,
                                           n_classes = soil_classes)
                      }),
                      filename  = file.path(strs, 'dsmartr_probabilities.tif'),
                      datatype  = 'FLT4S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)

  # think of ordered as the names of sorted
  # e.g. [1] 3 2 1 4 or [1]  3   2    1 4
  #          8 5 0 0         0.4 0.25 0 0

  message(paste0(Sys.time(),
                 ': Probabilities ordered.'))
  endCluster()

  # Assemble outputs
  tallied_predictions <- raster::stack(file.path(strs, 'tallied_predictions.tif'))
  # faster to go lookup$CLASS but then can't be flexible about colnames
  names(tallied_predictions) <- lookup[[2]]

  tallied_probabilities <-
    raster::stack(file.path(strs, 'tallied_probabilities.tif'))
  names(tallied_probabilities) <- lookup[[2]]

  dsmartr_predictions <-
    raster::stack(file.path(strs, 'dsmartr_predictions.tif'))
  # give these rasters categories for these values
  # nb seq(layers(x)) is not faster than 1:n in this case - 3x slower
  # no danger of trying to get length of missing object; fn would fail before
  # this point in the case of a missing file
  for (i in seq(nlayers(dsmartr_predictions))) {
    levels(dsmartr_predictions[[i]]) <- lookup
  }
  names(dsmartr_predictions) <- paste0('mostlikely_',
                                       seq(nlayers(dsmartr_predictions)))
  dsmartr_probabilities <-
    raster::stack(file.path(strs, 'dsmartr_probabilities.tif'))
  names(dsmartr_probabilities) <- paste0('mostlikely_prob_',
                                         seq(nlayers(dsmartr_probabilities)))

  collated_outputs <- list("dsmartr_predictions"   = dsmartr_predictions,
                           "dsmartr_probabilities" = dsmartr_probabilities,
                           "tallied_predictions"   = tallied_predictions,
                           "tallied_probabilities" = tallied_probabilities)
  rm(list = c('soil_classes', 'n_iterations'), envir = parent.frame())
  message(paste0(Sys.time(),
                 ': Collation complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'tallies')))
  collated_outputs
}
