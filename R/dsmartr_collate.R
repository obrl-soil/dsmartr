#' Collate dsmartr iterations
#'
#' Processes the outputs of \code{\link{dsmartr_iterate}}
#' @param iteration_stack RasterStack or Brick; output \code{iteration_maps} from
#' \code{\link{dsmartr_iterate}}.
#' @param lookup Data Frame; contains raster values and corresponding soil class labels. Example
#' \code{str(heronvale_soilnames)}:
#'  \preformatted{'data.frame':	36 obs. of  3 variables:
#' $ ID   : int  1 2 3 4 5 6 7 8 9 10 ...
#' $ CLASS: chr  "Ad" "An" "Bb" "Bh" ...
#' $ NAME : chr  "Andromache" "Andergrove" "Balberra" "Benholme" ...}
#' There's no need to have both CLASS and NAME fields but it can come in handy for map-making etc.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \itemize{
#' \item{\code{dsmartr_predictions}: RasterStack, soil class prediction maps in order of most to
#' least-probable}
#' \item{\code{dsmartr_probabilties}: RasterStack, a probability surface for each layer of
#' \code{dsmartr_predictions}.}
#' \item{\code{tallied_predictions}: Optional RasterStack containing per-pixel tallies of soil
#' class occurrence.}
#' \item{\code{tallied_probabilities}: Optional RasterStack containing per-pixel probabilities of
#' soil class occurrence}.}
#'  All outputs are written to disk as multiband GeoTIFFs before being assigned to the global
#'  environment.
#' @note This function can generate very large R temporary files. A rough minimum disk space
#' requirement is \code{n} iterations * \code{n} covariate cells * 8 bytes * 2.
#' @examples \dontrun{
#' # run dsmartr_iterate() with the example data, then:
#' LUT <- levels(iteration_maps[[1]])[[1]]
#' LUT <- LUT[!(LUT$ID == 0), ]
#' dsmartr_collate(iteration_stack = iteration_maps, lookup = LUT,
#' cpus = max(1, (parallel::detectCores() - 1)))}
#' @importFrom raster beginCluster calc clusterR endCluster nlayers stack
#' @importFrom stats na.omit
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
dsmartr_collate <- function(iteration_stack = NULL,
                            lookup          = NULL,
                            cpus            = 1) {

  if (!dir.exists("tallies")) {
    dir.create("tallies", showWarnings = F)
  }
  strs         <- file.path(getwd(), 'tallies')
  classes      <- nrow(lookup)
  n_iterations <- nlayers(iteration_stack)

  ### cell by cell calc functions ###

  # class_count produces an integer vector counting the number of times a given soil was
  # predicted on this pixel
  # e.g. 0 5 8 0 == soil 1 was not predicted, soil 2 was predicted 5 times, soil 3 was predicted
  # 8 times, soil 4 was not predicted.
  class_count <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else {
      tabulate(x, nbins = classes)
    }
  }

  # class_prop takes the vector produced by class_count and normalises it against the
  # total number of model runs
  # e.g. for 20 runs, 0 5 8 0 becomes 0 0.25 0.4 0
  class_prop  <- function(x) {round((x / n_iterations), 3)}

  # stack_cell_sort sorts the elements of class_order or probs_order from largest to smallest.
  # e.g. 0 0.25 0.4 0 becomes [1] 0.4 0.25 0 0
  stack_cell_sort <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else {
      sort(x, decreasing = TRUE, na.last = TRUE)
    }
  }

  # stack_cell_order returns a vector of the position of the elements of probs_order
  # from largest to smallest.
  # e.g. 0 0.25 0.4 0 becomes [1] 3 2 1 4
  # these values correspond to the lookup ID column, so they can be linked to soil class name
  stack_cell_order <- function(x) {
    if (is.na(sum(x))) {
      rep(NA, classes)
    } else {
      order(x, decreasing = TRUE, na.last = TRUE)
    }
  }

  # note that its easiest to think about the results of stack_cell_order as the names
  # of stack_cell_sort  e.g. [1] 3 2 1 4 or [1]  3   2   1 4
  #                              8 5 0 0        0.4 0.25 0 0

  ###
  beginCluster(cpus)
  assign('classes',    classes,    envir = .GlobalEnv)
  assign('iterations', n_iterations, envir = .GlobalEnv)

  message(paste0(Sys.time(), ': dsmartr collation in progress...'))

  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)

  # 1. get a stack of predicted class counts
  counts <- clusterR(x         = iteration_stack,
                     fun       = calc,
                     args      = list(fun = class_count),
                     export    = 'classes',
                     filename  = file.path(strs, 'tallied_predictions.tif'),
                     datatype  = 'INT2S',
                     format    = 'GTiff',
                     NAflag    = -9999,
                     overwrite = TRUE)

  counts <- raster::stack(file.path(strs, 'tallied_predictions.tif'))
  names(counts) <- as.vector(lookup[ , 2])
  assign('tallied_predictions', counts, envir = .GlobalEnv)

  setTxtProgressBar(pb, 25)

  # 2. express that stack as a probability
  probs <- clusterR(x         = counts,
                    fun       = calc,
                    args      = list(fun = class_prop),
                    export    = 'iterations',
                    filename  = file.path(strs, 'tallied_probabilities.tif'),
                    datatype  = 'FLT4S',
                    format    = 'GTiff',
                    NAflag    = -9999,
                    overwrite = TRUE)
  probs <- raster::stack(file.path(strs, 'tallied_probabilities.tif'))
  names(probs) <- as.vector(lookup[ , 2])
  assign('tallied_probabilities', probs, envir = .GlobalEnv)

  setTxtProgressBar(pb, 50)

  # 3. order probs by most to least probable
  ordered <- clusterR(x         = na.omit(probs),
                      fun       = calc,
                      args      = list(fun = stack_cell_order),
                      export    = 'classes',
                      filename  = file.path(strs, 'dsmart_predictions.tif'),
                      datatype  = 'INT2S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)

  ordered <- raster::stack(file.path(strs, 'dsmart_predictions.tif'))
  # give these rasters categories for these values
  for (i in 1:nlayers(ordered)) {
    levels(ordered[[i]]) <- lookup
  }
  names(ordered) <- paste0('mostlikely_', 1:nlayers(ordered))
  assign('dsmart_predictions', ordered, envir = .GlobalEnv)

  setTxtProgressBar(pb, 75)

  # 4. sort probs by most to least probable
  sorted  <- clusterR(x         = probs,
                      fun       = calc,
                      args      = list(fun = stack_cell_sort),
                      export    = 'classes',
                      filename  = file.path(strs, 'dsmart_probabilities.tif'),
                      datatype  = 'FLT4S',
                      format    = 'GTiff',
                      NAflag    = -9999,
                      overwrite = TRUE)
  sorted <- raster::stack(file.path(strs, 'dsmart_probabilities.tif'))
  names(sorted) <- paste0('mostlikely_prob_', 1:nlayers(sorted))
  assign('dsmart_probabilities', sorted, envir = .GlobalEnv)

  setTxtProgressBar(pb, 100)
  close(pb)
  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'tallies')))
  endCluster()
  rm(iterations, pos = .GlobalEnv)
  rm(classes, pos = .GlobalEnv)
}
