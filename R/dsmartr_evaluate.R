#' Calculate probability gap
#'
#' Calculates and maps the difference between the first and second most-probable
#' dsmartr probability surfaces. Requires outputs of \code{\link{collate}}.
#' @param dsmartr_probs RasterBrick; 'dsmartr_probabilities' output by
#'   \code{\link{collate}}. Alternatively, probability maps output by
#'   \code{\link{most_likely}} can be used, or a list of two rasters read from disk.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{probability_gap}: RasterLayer depicting the probability gap. Written to disk as
#'   GeoTIFF.
#' @note This function is often called the 'confusion index', but has been renamed as that term is
#'   used in multiple contexts within the scientific literature.
#' @examples \dontrun{
#' # run collate() with the example data then:
#' pgap1 <- eval_pgap(dsmartr_probs = collated[['dsmartr_probabilities']][[1:2]],
#' cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # or supply unstacked maps after running unstack() (slightly faster)
#' pgap2 <- eval_pgap(dsmartr_probs = most_likely_soil[c('most_likely_prob_1',
#'                                                       'most_likely_prob_2')],
#'                    cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # or read from file
#' mpp_1 <- raster(file.path(getwd(), 'most_likely_maps', 'most_likely_1.tif'))
#' mpp_2 <- raster(file.path(getwd(), 'most_likely_maps', 'most_likely_2.tif'))
#' pgap3 <- eval_pgap(dsmartr_probs = list(mpp_1, mpp_2),
#'                    cpus = max(1, (parallel::detectCores() - 1)))
#' }
#' @importFrom methods is
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
eval_pgap <- function(dsmartr_probs = NULL, cpus = 1) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  message(paste0(Sys.time(), ': dsmartr probability gap calculation in progress...'))
  dsmartr_probs <- if(is(dsmartr_probs, 'list')) {
    raster::stack(dsmartr_probs)
  } else { dsmartr_probs }

  prob_gap <- function(cell = NULL) {
    1 - (cell[[1]] - cell[[2]])
    }
  beginCluster(cpus)
  probability_gap <- clusterR(x         = dsmartr_probs,
                              fun       = prob_gap,
                              filename  = file.path(strs, 'probability_gap.tif'),
                              datatype  = 'FLT4S',
                              NAflag    = -9999,
                              overwrite = TRUE)
  endCluster()
  message(paste0(Sys.time(), ': ...complete. Function outputs can be located at ',
                 file.path(getwd(), 'evaluation')))
  probability_gap
}

#' Number of classes predicted
#'
#' Counts how many different soil classes were predicted on a pixel
#' @keywords internal
#' @param input A vector of integers n soil classes long, where the integers refer to the number of
#' times a soil class was predicted.
#' @param n_iterations Integer; the total number of model runs.
#' @param noise_cutoff Decimal; proportion of predictions to be considered 'noise' and ignored.
#' Allowable values range between 0 and 1.
#' @return Integer; the number of distinct soil classes predicted on the pixel (optionally, above
#'  a supplied threshold).
#' @note This is a helper function and not widely applicable. It is expected to be used in concert
#' with \code{\link[raster:calc]{calc}}, applied to a stack of tallied soil class predictions
#' produced by \code{\link[dsmartr:eval_npred]{eval_npred}}.
#' @examples \dontrun{
#' # a set of tallies from the heronvale demo data
#' counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
#' 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0))
#' n_soils <- n_soils_predicted(counts_19291)
#' n_soils_10 <- n_soils_predicted(counts_19291, 100, 0.1)
#' }
#'
n_predicted <- function(input = NULL, n_iterations = NULL, noise_cutoff = NULL) {

  if(!is.null(noise_cutoff)) {
    stopifnot(in_range(noise_cutoff, 0, 1, strict = FALSE))
  }

  tol <- if(!is.null(noise_cutoff)) {
    ceiling(n_iterations * noise_cutoff)
  } else {
    0L
  }

  if(all(is.na(input))) {
    NA_integer_
    } else {
      length(input[input > tol]) # NB NA > 0 but there should never be an NA fed to this fn
    }
}

#' Calculate soil classes per pixel
#'
#' Calculates and maps the number of distinct soil classes predicted. Requires outputs of
#' \code{\link[dsmartr:collate]{collate}}.
#' @param tallied_preds RasterBrick; 'tallied_predictions' output by
#' \code{\link[dsmartr:collate]{collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @param n_iterations Optional Integer; the number of iterations supplied to
#' \code{\link[dsmartr:iterate]{iterate}}.
#' @param noise_cutoff Optional Decimal; proportion of predictions to be considered 'noise' and
#'  ignored. Acceptable values range between 0 and 1 inclusive. Defaults to 0.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of distinct soils predicted
#'   per pixel more than \code{n_iterations * noise_cutoff} times. Written to disk as GeoTIFF.
#' @note Fewer classes predicted on a pixel generally indicates higher internal model confidence at
#'   that location.
#' @examples \dontrun{
#' # run iterate() and collate() with the example data then:
#' nxpred <- eval_npred(tallied_preds = collated[['tallied_predictions']],
#'                      cpus = max(1, (parallel::detectCores() - 1)),
#'                      n_iterations = nlayers(iteration_maps), noise_cutoff = 0.1)
#' }
#' @import parallel
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
eval_npred <- function(tallied_preds = NULL,
                        cpus          = 1,
                        n_iterations  = NULL,
                        noise_cutoff  = NULL) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  message(paste0(Sys.time(), ': dsmartr notable soils predicted calculation in progress...'))
  beginCluster(cpus)
  assign('noise_cutoff', noise_cutoff, envir = parent.frame())
  assign('n_iterations', n_iterations, envir = parent.frame())

  n_classes_predicted_x <- clusterR(tallied_preds,
                                    fun       = calc,
                                    args = list(fun = function(cell) {
                                      n_predicted(cell, n_iterations, noise_cutoff)
                                    }),
                                    filename  = file.path(strs,
                                                          paste0('n_classes_predicted_over_',
                                                                 (n_iterations * noise_cutoff),
                                                                 '.tif')),
                                    datatype  = 'INT2S',
                                    NAflag    = -9999,
                                    overwrite = TRUE)

  endCluster()
  rm(list = c('n_iterations', 'noise_cutoff'), envir = parent.frame())
  message(paste0(Sys.time(), ': ...complete. Function outputs can be located at ',
                 file.path(getwd(), 'evaluation')))
  n_classes_predicted_x
}

#' Number of tied soil classes
#'
#' Checks existence of ties for most-probable soil and reports how many classes are tied on a given
#' pixel.
#' @keywords internal
#' @param input A vector of integers n soil classes long, where the integers refer to the number of
#' times a soil class was predicted.
#' @param n_iterations Integer; the total number of model runs.
#' @param noise_cutoff Decimal; proportion of predictions to be considered 'noise' and ignored.
#' Allowable values range between 0 and 1.
#' @return Integer; the number of distinct soil classes predicted on the pixel (optionally, above
#'  a supplied threshold).
#' @note This is a helper function and not widely applicable. It is expected to be used in concert
#' with \code{\link[raster:calc]{calc}}, applied to a stack of tallied soil class predictions
#' produced by \code{\link[dsmartr:eval_npred]{eval_npred}}.
#' @examples \dontrun{
#' # a set of tallies from the heronvale demo data
#' counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
#' 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0))
#' ties <- tie_finder(counts_19291) # ties exist, but not for most-probable
#' }
#'
tie_finder <- function(input = NULL) {
  if (all(is.na(input))) {
    NA_integer_
  } else {
    y <- max(input, na.rm = TRUE)
    z <- length(input[input == y])
    if(z == 1L) { 0L } else { z }
  }
}

#' Detect ties for most-probable soil
#'
#' Calculates and maps out locations where the most-probable soil is a tie between two or more
#' classes. Requires output from \code{\link{collate}}.
#' @param tallied_preds RasterBrick; \code{tallied_predictions} output by
#'   \code{\link{collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of soils tied for
#'   most-probable per pixel. Written to disk as GeoTIFF.
#' @examples \dontrun{
#' # run iterate() and collate() with the example data then:
#' tiemap <- eval_ties(tallied_preds = collated[['tallied_predictions']],
#'                    cpus = max(1, (parallel::detectCores() - 1)))
#'                    }
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
eval_ties <- function(tallied_preds = NULL,
                      cpus          = 1) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  message(paste0(Sys.time(), ': dsmartr tie-finder calculation in progress...'))
  beginCluster(cpus)
  tie_map <- clusterR(tallied_preds,
                      fun = calc,
                      args = list(fun = function(cell) {
                        tie_finder(input = cell)
                        }),
                      filename  = file.path(strs, 'most_probable_ties.tif'),
                      datatype  = 'INT2S',
                      NAflag    = -9999,
                      overwrite = TRUE)

  endCluster()
  message(paste0(Sys.time(), ': ...complete. Function outputs can be located at ',
                 file.path(getwd(), 'evaluation')))
  tie_map
}
