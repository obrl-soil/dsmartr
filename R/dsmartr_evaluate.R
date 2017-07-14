#' Calculate probability gap
#'
#' Calculates and maps the difference between the first and second most-probable
#' dsmartr probability surfaces. Requires outputs of \code{\link{dsmartr_collate}}.
#' @param dsmartr_probs RasterBrick; 'dsmartr_probabilities' output by
#'  \code{\link{dsmartr_collate}}. Alternatively, probability maps output by
#'  \code{\link{dsmartr_unstack}} can be used, or a list of two rasters read from disk.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{probability_gap}: RasterLayer depicting the probability gap. Written to disk as
#' GeoTIFF before being assigned to the global environment.
#' @note This function is often called the 'confusion index', but has been renamed as that term is
#' used in multiple contexts within the scientific literature.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data then:
#' dsmartr_eval_pgap(dsmartr_probs = dsmartr_probabilities, max(1, (parallel::detectCores() - 1)))
#'
#' # or supply unstacked maps
#' dsmartr_eval_pgap(dsmartr_probs = most_prob_prob_maps[1:2],
#' max(1, (parallel::detectCores() - 1)))
#'
#' # or read from file
#' mpp_1 <- raster(file.path(getwd(), 'most_probable_maps', 'mostlikely_1.tif'))
#' mpp_2 <- raster(file.path(getwd(), 'most_probable_maps', 'mostlikely_2.tif'))
#' dsmartr_eval_pgap(dsmartr_probs = list(mpp_1, mpp_2),
#' max(1, (parallel::detectCores() - 1)))
#' }
#' @export
dsmartr_eval_pgap <- function(dsmartr_probs = NULL, cpus = 1) {

  if (!dir.exists('evaluation')) {
    dir.create('evaluation', showWarnings = F)
  }
  strs <- file.path(getwd(), 'evaluation')

  dsmartr_probs <- if(is(dsmartr_probs, 'list')) {
    raster::stack(dsmartr_probs)
  } else { dsmartr_probs }

  prob_gap  <- function(x) {1 - (x[[1]] - x[[2]])}
  beginCluster(cpus)
  pgap  <- clusterR(na.omit(dsmartr_probs),
                    fun       = prob_gap,
                    filename  = file.path(strs, 'probability_gap.tif'),
                    datatype  = 'FLT4S',
                    NAflag    = -9999,
                    overwrite = TRUE)
  assign('probability_gap', pgap, envir = .GlobalEnv)
  endCluster()

}

#' Calculate soil classes per pixel
#'
#' Calculates and maps the number of distinct soil classes predicted by dsmartr across \code{n}
#' iterations. Requires outputs of \code{\link{dsmartr_collate}} (with option
#' \code{keep_tallies = TRUE}).
#' @param tallied_preds RasterBrick; \code{tallied_predictions} output by
#' \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @param n_iterations Integer; the number of iterations that weresupplied to
#' \code{\link{dsmartr_iterate}}.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of distinct soils predicted
#' per pixel. Written to disk as GeoTIFF before being assigned to the global environment.
#' @note Fewer classes predicted on a pixel generally indicates higher internal model confidence at
#' that location.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data then:
#' dsmartr_eval_npred(tallied_preds = tallied_predictions, max(1, (parallel::detectCores() - 1)),
#' n_iterations = nlayers(iteration_maps))}
#' @export
dsmartr_eval_npred <- function(tallied_preds = NULL,
                               cpus          = 1,
                               n_iterations  = NULL) {

  if (is.null(n_iterations)) {
    stop('total number of DSMART model iterations must be supplied.')
  }

  if (!dir.exists('evaluation')) {
    dir.create('evaluation', showWarnings = F)
  }
  strs <- file.path(getwd(), 'evaluation')

  ### per-cell functions

  n_classes_predicted <- function(x) {
    ifelse(is.na(sum(x)), NA, length(x[x > 0]))
  }

  beginCluster(cpus)
  assign('n_iterations', n_iterations, envir = .GlobalEnv)

  npred <- clusterR(tallied_preds,
                    fun       = calc,
                    args      = list(fun = n_classes_predicted),
                    filename  = file.path(strs, 'n_classes_predicted.tif'),
                    export    = 'n_iterations',
                    NAflag    = -9999,
                    datatype  = 'INT2S',
                    overwrite = TRUE)

  assign('n_classes_predicted', npred, envir = .GlobalEnv)
  endCluster()
  rm(n_iterations)

}

#' Calculate notable soil classes per pixel
#'
#' Calculates and maps the number of distinct soil classes predicted in more than a
#' given proportion of iterations. Requires outputs of \code{\link{dsmartr_collate()}} (with option
#' \code{keep_tallies = TRUE}).
#' @param tallied_preds RasterBrick; 'tallied_predictions' output by \code{\link{dsmartr_collate()}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @param n_iterations Integer; the number of iterations supplied to \code{\link{dsmartr_iterate()}}.
#' @param noise_cutoff Decimal; proportion of predictions to be considered 'noise' and ignored.
#' @return \code{n_classes_predicted}: RasterLayer depicting the number of distinct soils predicted
#' per pixel more than \code{n_iterations * noise_cutoff} times. Written to disk as GeoTIFF before
#' being assigned to the global environment.
#' @note Fewer classes predicted on a pixel generally indicates higher internal model confidence at
#' that location.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data then:
#' dsmartr_eval_nxpred(tallied_preds = tallied_predictions, max(1, (parallel::detectCores() - 1)),
#' n_iterations = nlayers(iteration_maps), noise_cutoff = 0.1)}
#' @export
dsmartr_eval_nxpred <- function(tallied_preds = NULL,
                               cpus          = 1,
                               n_iterations  = NULL,
                               noise_cutoff  = 0.1) {

  if (is.null(n_iterations)) {
    stop('total number of DSMART model iterations must be supplied.')
  }

  if (!dir.exists('evaluation')) {
    dir.create('evaluation', showWarnings = F)
  }
  strs <- file.path(getwd(), 'evaluation')

  n_classes_predicted_x <- function(x) {
    ifelse(is.na(sum(x)), NA, length(x[x > (ceiling(n_iterations * noise_cutoff))]))
  }

  beginCluster(cpus)
  assign('noise_cutoff', noise_cutoff, envir = .GlobalEnv)
  assign('n_iterations', n_iterations, envir = .GlobalEnv)

  npredx <- clusterR(tallied_preds,
                     fun       = calc, args = list(fun = n_classes_predicted_x),
                     export    = c('n_iterations', 'noise_cutoff'),
                     filename  = file.path(strs,
                                           paste0('n_classes_predicted_over_',
                                                  (n_iterations * noise_cutoff),
                                                  '.tif')),
                     datatype  = 'INT2S',
                     NAflag    = -9999,
                     overwrite = TRUE)

  assign('n_classes_predicted_x', npredx, envir = .GlobalEnv)
  endCluster()
  rm(noise_cutoff)
  rm(n_iterations)

}
