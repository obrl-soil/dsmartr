#' Collate dsmartr iterations
#'
#' Processes the outputs of \code{\link{dsmartr_iterate}}
#' @param iteration_stack RasterStack or Brick; output \code{iteration_maps} from
#'   \code{\link{dsmartr_iterate}}.
#' @param lookup Data Frame; contains raster values and corresponding soil class labels. Example:
#'   \preformatted{'data.frame':	38 obs. of  2 variables:
#'   $ VALUE: int  1 2 3 4 5 6 7 8 9 10 ...
#'   $ CLASS: chr  "Ad" "An" "Bb" "Bh" ...}
#' @param cpus Integer; number of processors to use in parallel.
#' @return A list of four RasterStacks:
#' \itemize{
#'   \item{\code{dsmartr_predictions}: RasterStack, soil class prediction maps in order of most to
#'   least-probable}
#'   \item{\code{dsmartr_probabilties}: RasterStack, a probability surface for each layer of
#'   \code{dsmartr_predictions}.}
#'   \item{\code{tallied_predictions}: Optional RasterStack containing per-pixel tallies of soil
#'   class occurrence.}
#'   \item{\code{tallied_probabilities}: Optional RasterStack containing per-pixel probabilities of
#'   soil class occurrence}.}
#'
#' All outputs are written to disk as multiband GeoTIFFs before being assigned to the global
#' environment.
#' @note This function can generate very large R temporary files. A rough minimum disk space
#'   requirement is \code{n} iterations * \code{n} covariate cells * 8 bytes * 2.
#' @examples \dontrun{
#' # run dsmartr_iterate() with the example data, then:
#' LUT <- levels(iteration_maps[[1]])[[1]]
#' collated <- dsmartr_collate(iteration_stack = iteration_maps, lookup = LUT,
#' cpus = max(1, (parallel::detectCores() - 1)))}
#' @importFrom raster beginCluster calc clusterR endCluster getCluster nlayers stack
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
#'
dsmartr_collate <- function(iteration_stack = NULL,
                            lookup          = NULL,
                            cpus            = 1) {
  message(paste0(Sys.time(), ': dsmartr collation in progress...'))

  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  if (!dir.exists("tallies")) {
    dir.create("tallies", showWarnings = F)
  }
  strs           <- file.path(getwd(), 'tallies')
  soil_classes   <- nrow(lookup)
  n_iterations   <- nlayers(iteration_stack)

  beginCluster(n = cpus)
  assign('soil_classes', soil_classes, envir = parent.frame())
  assign('n_iterations', n_iterations, envir = parent.frame())
  # 1. tally up how many times each soil class was predicted
  counts <- clusterR(x         = iteration_stack,
                     fun       = calc,
                     args      = list(fun =
                                        # produces an integer vector counting the number of times a
                                        # given soil was predicted on this pixel
                                        # e.g. 0 5 8 0 == soil 1 was not predicted, soil 2 was
                                        # predicted 5 times, soil 3 was predicted 8 times, soil 4
                                        # was not predicted.
                                        function(x) {
                                          if (is.na(sum(x))) {
                                            rep(NA, soil_classes)
                                          } else {
                                            tabulate(x, nbins = soil_classes)
                                          }
                                        }),
                     export    = c('soil_classes'),
                     filename  = file.path(strs, 'tallied_predictions.tif'),
                     datatype  = 'INT2S',
                     format    = 'GTiff',
                     NAflag    = -9999,
                     overwrite = TRUE)

  setTxtProgressBar(pb, 30)

  # 2. express the tallies as probablities
  probs <- clusterR(x         = counts,
                    fun       = calc,
                    args      = list(fun =
                                     # class_prop takes the vector produced by class_count and
                                     # normalises it against the total number of model runs
                                     # e.g. for 20 runs, 0 5 8 0 becomes 0 0.25 0.4 0
                                     function(x) {round((x / n_iterations), 3)}),
                    export    = c('n_iterations'),
                    filename  = file.path(strs, 'tallied_probabilities.tif'),
                    datatype  = 'FLT4S',
                    format    = 'GTiff',
                    NAflag    = -9999,
                    overwrite = TRUE)

  setTxtProgressBar(pb, 60)

  # 3. get the order of counts
  ordered <- clusterR(x         = counts,
                      fun       = calc,
                      args      = list(fun =
                                       # returns a vector of the position of the elements of
                                       # counts from largest to smallest.
                                       # e.g. 0 5 8 0 becomes [1] 3 2 1 4
                                       # these values correspond to the lookup ID column, so they
                                       # can be linked to soil class name
                                       # where ties for x-most-probable exist, they are shuffled
                                       # rather than being output in ascending order within the tie
                                  function(x) {
                                         if (is.na(sum(x))) {
                                           rep(NA, soil_classes)
                                         } else {
                                             order(x, sample(x, size = length(x), replace = FALSE),
                                                   decreasing = TRUE, na.last = TRUE)
                                           }}),
                      export    = c('soil_classes'),
                      filename  = file.path(strs, 'dsmartr_predictions.tif'),
                      datatype  = 'INT2S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)

  setTxtProgressBar(pb, 60)

  # 4. sort probs by most to least probable
  sorted  <- clusterR(x         = probs,
                      fun       = calc,
                      args      = list(fun =
                                       # sorts the elements of class_order or probs_order from
                                       # largest to smallest.
                                       # e.g. 0 0.25 0.4 0 becomes [1] 0.4 0.25 0 0
                                       # still works with the tie-shuffle above as shuffled soils
                                       # are equally probable
                                       function(x) {
                                         if (is.na(sum(x))) {
                                           rep(NA, soil_classes)
                                         } else {
                                           sort(x, decreasing = TRUE, na.last = TRUE)
                                         }
                                       }),
                      export    = c('soil_classes'),
                      filename  = file.path(strs, 'dsmartr_probabilities.tif'),
                      datatype  = 'FLT4S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)

  # think of ordered as the names of sorted
  # e.g. [1] 3 2 1 4 or [1]  3   2   1 4
  #         8 5 0 0        0.4 0.25 0 0

  setTxtProgressBar(pb, 90)
  endCluster()

  # Assemble outputs
  tallied_predictions <- raster::stack(file.path(strs, 'tallied_predictions.tif'))
  names(tallied_predictions) <- as.vector(lookup[ , 2])

  tallied_probabilities <- raster::stack(file.path(strs, 'tallied_probabilities.tif'))
  names(tallied_probabilities) <- as.vector(lookup[ , 2])

  dsmartr_predictions <- raster::stack(file.path(strs, 'dsmartr_predictions.tif'))
  # give these rasters categories for these values
  for (i in 1:nlayers(dsmartr_predictions)) {
    levels(dsmartr_predictions[[i]]) <- lookup
  }
  names(dsmartr_predictions) <- paste0('mostlikely_', 1:nlayers(dsmartr_predictions))

  dsmartr_probabilities <- raster::stack(file.path(strs, 'dsmartr_probabilities.tif'))
  names(dsmartr_probabilities) <- paste0('mostlikely_prob_', 1:nlayers(dsmartr_probabilities))

  collated_outputs <- list("dsmartr_predictions"   = dsmartr_predictions,
                           "dsmartr_probabilities" = dsmartr_probabilities,
                           "tallied_predictions"   = tallied_predictions,
                           "tallied_probabilities" = tallied_probabilities)
  setTxtProgressBar(pb, 100)
  close(pb)
  rm(list = c('soil_classes', 'n_iterations'), envir = parent.frame())
  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'tallies')))
  collated_outputs
}
