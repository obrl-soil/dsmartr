#' Calculate probability gap
#'
#' Calculates and maps the difference between the first and second most-probable
#' dsmartr probability surfaces. Requires outputs of \code{\link{dsmartr_collate}}.
#' @param dsmartr_probs RasterBrick; 'dsmartr_probabilities' output by
#'   \code{\link{dsmartr_collate}}. Alternatively, probability maps output by
#'   \code{\link{dsmartr_most_likely}} can be used, or a list of two rasters read from disk.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{probability_gap}: RasterLayer depicting the probability gap. Written to disk as
#'   GeoTIFF.
#' @note This function is often called the 'confusion index', but has been renamed as that term is
#'   used in multiple contexts within the scientific literature.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data then:
#' pgap1 <- dsmartr_eval_pgap(dsmartr_probs = collated[['dsmartr_probabilities']][1:2],
#' cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # or supply unstacked maps (slightly faster)
#' pgap2 <- dsmartr_eval_pgap(dsmartr_probs = most_likely[c('most_likely_prob_1', 'most_likely_prob_2')],
#' #' cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # or read from file
#' mpp_1 <- raster(file.path(getwd(), 'most_likely_maps', 'most_likely_1.tif'))
#' mpp_2 <- raster(file.path(getwd(), 'most_likely_maps', 'most_likely_2.tif'))
#' pgap3 <- dsmartr_eval_pgap(dsmartr_probs = list(mpp_1, mpp_2),
#' cpus = max(1, (parallel::detectCores() - 1)))
#' }
#' @importFrom methods is
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
dsmartr_eval_pgap <- function(dsmartr_probs = NULL, cpus = 1) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  message(paste0(Sys.time(), ': dsmartr probability gap calculation in progress...'))
  dsmartr_probs <- if(is(dsmartr_probs, 'list')) {
    raster::stack(dsmartr_probs)
  } else { dsmartr_probs }

  prob_gap  <- function(x) {1 - (x[[1]] - x[[2]])}
  beginCluster(cpus)
  probability_gap <- clusterR(na.omit(dsmartr_probs),
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

#' Calculate soil classes per pixel
#'
#' Calculates and maps the number of distinct soil classes predicted by dsmartr across \code{n}
#' iterations. Requires output from \code{\link{dsmartr_collate}}.
#' @param tallied_preds RasterBrick; \code{tallied_predictions} output by
#'   \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of distinct soils predicted
#'   per pixel. Written to disk as GeoTIFF.
#' @note Fewer classes predicted on a pixel generally indicates higher internal model confidence at
#'   that location.
#' @examples \dontrun{
#' # run dsmartr_iterate() and dsmartr_collate() with the example data then:
#' npred <- dsmartr_eval_npred(tallied_preds = collated[['tallied_predictions']],
#' cpus = max(1, (parallel::detectCores() - 1)))}
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
dsmartr_eval_npred <- function(tallied_preds = NULL,
                               cpus          = 1) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')
  message(paste0(Sys.time(), ': dsmartr n soils predicted calculation in progress...'))
  beginCluster(cpus)
  n_classes_predicted <- clusterR(tallied_preds,
                                  fun       = calc,
                                  args      = list(fun = function(x) {
                                    ifelse(is.na(sum(x)), NA, length(x[x > 0]))
                                  }),
                                  filename  = file.path(strs, 'n_classes_predicted.tif'),
                                  NAflag    = -9999,
                                  datatype  = 'INT2S',
                                  overwrite = TRUE)

  endCluster()
  message(paste0(Sys.time(), ': ...complete. Function outputs can be located at ',
                 file.path(getwd(), 'evaluation')))
  n_classes_predicted

}

#' Calculate notable soil classes per pixel
#'
#' Calculates and maps the number of distinct soil classes predicted in more than a
#' given proportion of iterations. Requires outputs of \code{\link{dsmartr_collate}} (with option
#' \code{keep_tallies = TRUE}).
#' @param tallied_preds RasterBrick; 'tallied_predictions' output by \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @param n_iterations Integer; the number of iterations supplied to \code{\link{dsmartr_iterate}}.
#' @param noise_cutoff Decimal; proportion of predictions to be considered 'noise' and ignored.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of distinct soils predicted
#'   per pixel more than \code{n_iterations * noise_cutoff} times. Written to disk as GeoTIFF.
#' @note Fewer classes predicted on a pixel generally indicates higher internal model confidence at
#'   that location.
#' @examples \dontrun{
#' # run dsmartr_iterate() and dsmartr_collate() with the example data then:
#' nxpred <- dsmartr_eval_nxpred(tallied_preds = collated[['tallied_predictions']],
#' cpus = max(1, (parallel::detectCores() - 1)),
#' n_iterations = nlayers(iteration_maps), noise_cutoff = 0.1)}
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
dsmartr_eval_nxpred <- function(tallied_preds = NULL,
                                cpus          = 1,
                                n_iterations  = NULL,
                                noise_cutoff  = 0.1) {

  if (is.null(n_iterations)) {
    stop('total number of DSMART model iterations must be supplied.')
  }

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  n_classes_predicted_x <- function(x) {
    ifelse(is.na(sum(x)), NA, length(x[x > (ceiling(n_iterations * noise_cutoff))]))
  }

  message(paste0(Sys.time(), ': dsmartr notable soils predicted calculation in progress...'))
  beginCluster(cpus)
  assign('noise_cutoff', noise_cutoff, envir = parent.frame())
  assign('n_iterations', n_iterations, envir = parent.frame())

  n_classes_predicted_x <- clusterR(tallied_preds,
                                    fun       = calc, args = list(fun = n_classes_predicted_x),
                                    export    = c('n_iterations', 'noise_cutoff'),
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


#' Detect ties for most-probable soil
#'
#' Calculates and maps out locations where the most-probable soil is a tie between two or more
#' classes. Requires output from \code{\link{dsmartr_collate}}.
#' @param tallied_preds RasterBrick; \code{tallied_predictions} output by
#'   \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of soils tied for
#'   most-probable per pixel. Written to disk as GeoTIFF.
#' @examples \dontrun{
#' # run dsmartr_iterate() and dsmartr_collate() with the example data then:
#' tiemap <- dsmartr_eval_ties(tallied_preds = collated[['tallied_predictions']],
#' cpus = max(1, (parallel::detectCores() - 1)))}
#' @importFrom raster beginCluster calc clusterR endCluster writeRaster
#' @export
dsmartr_eval_ties <- function(tallied_preds = NULL,
                              cpus          = 1) {

  if (!dir.exists(file.path(getwd(), 'evaluation'))) {
    dir.create(file.path(getwd(), 'evaluation'), showWarnings = FALSE)
  }
  strs <- file.path(getwd(), 'evaluation')

  tie_finder <- function(x) {
    if (is.na(sum(x))) {
      NA
    } else {
      y <- max(x, na.rm = TRUE)
      z <- length(x[x==y])
      ifelse(z == 1, 0, z)
    }}

  message(paste0(Sys.time(), ': dsmartr tie-finder calculation in progress...'))
  beginCluster(cpus)
  tie_map <- clusterR(tallied_preds,
                      fun = calc,
                      args = list(fun = tie_finder),
                      filename  = file.path(strs, 'most_probable_ties.tif'),
                      datatype  = 'INT2S',
                      NAflag    = -9999,
                      overwrite = TRUE)

  endCluster()
  message(paste0(Sys.time(), ': ...complete. Function outputs can be located at ',
                 file.path(getwd(), 'evaluation')))
  tie_map
}
