#' Extract most-probable soil maps
#'
#' This function takes the outputs of \code{\link{dsmartr_collate}} and extracts the top \code{n}
#' most-probable soil maps. Optionally, their respective probability surfaces are also extracted.
#' @param dsmart_preds RasterBrick; output \code{dsmart_predictions} from
#'  \code{\link{dsmartr_collate}}.
#' @param dsmart_probs RasterBrick; output \code{dsmart_probabilities} from \code{\link{dsmartr_collate}};
#'  optional.
#' @param n_maps integer; the number of most-probable maps desired.
#' @return Two lists of RasterLayer objects:
#' \itemize{
#'   \item{\code{most_prob_maps}: A list of \code{n_maps} most-probable maps, from most to least
#'   probable.}
#'   \item{\code{most_prob_prob_maps}: Optional; A list of \code{n_maps} probability surfaces
#'   associated with \code{most_prob_maps}.}}
#'  All outputs are written to disk as GeoTIFFs before being assigned to the global environment.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data and then:
#' dsmartr_most_probable(dsmart_preds = dsmart_predictions, dsmart_probs = dsmart_probabilities,
#' n_maps = 2)}
#' @importFrom raster unstack ratify writeRaster
#' @export
dsmartr_most_probable <- function(dsmart_preds = NULL,
                                  dsmart_probs = NULL,
                                  n_maps       = 2) {
  if (!dir.exists('most_probable_maps/')) {
    dir.create('most_probable_maps/', showWarnings = F)
  }
  mp_dir <- file.path(getwd(), 'most_probable_maps')

  if (!is.null(dsmart_preds)) {
    suppressWarnings(map_list <- raster::unstack(dsmart_preds[[1:n_maps]]))
  } else {
    stop('dsmart_preds raster stack has not been specified')
  }

  lapply(map_list, function(i) ratify(i))

  outs_list <- mapply(FUN = function(x, i) {
    writeRaster(x,
                filename  = file.path(mp_dir, paste0('mostlikely_', i ,'.tif')),
                format    = 'GTiff',
                NAflag    = -9999,
                datatype  = 'INT2S',
                overwrite = TRUE)
  },
  x = map_list, i = seq_along(1:n_maps)
  )

  assign('most_prob_maps', outs_list, envir = .GlobalEnv)
  message(paste0(Sys.time(), ': ', n_maps, ' most-likely soils maps produced.'))

  ### optional probability surfaces
  if(!is.null(dsmart_probs)) {
    suppressWarnings(probmap_list <- raster::unstack(dsmart_probs[[1:n_maps]]))
  }

  probsouts_list <- mapply(FUN = function(x, i) {
    writeRaster(x,
                filename  = file.path(mp_dir, paste0('mostlikely_prob_', i, '.tif')),
                format    = 'GTiff',
                datatype  = 'FLT4S',
                NAflag    = -9999,
                overwrite = TRUE)
  },
  x = probmap_list,
  i = seq_along(1:n_maps)
  )

  assign('most_prob_prob_maps', probsouts_list, envir = .GlobalEnv)
  message(paste0(Sys.time(), ': ', n_maps, ' probability maps produced.'))

}

#' Produce soil class probability surfaces
#' @param tallied_probs RasterBrick; \code{tallied_probabilities} output by
#' \code{\link{dsmartr_collate}}.
#' @param soil_class String; Soil class(es) of interest, if only particular maps are desired.
#' @param lookup Data Frame; contains raster values and corresponding soil class labels. See
#' \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{class_maps}: List of RasterLayers; probability surfaces for each soil class.
#'
#' All outputs are written to disk as GeoTIFFs before being assigned to the global environment.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data and then
#'
#' # all classes:
#' dsmartr_class_maps(tallied_probs = tallied_probabilities, lookup = LUT)
#'
#' # just map two classes of interest:
#' dsmartr_class_maps(tallied_probs = tallied_probabilities,
#' soil_class = c('BL', 'CO'), lookup = LUT)}
#' @importFrom raster unstack writeRaster
#' @export
dsmartr_class_maps <- function(tallied_probs = NULL,
                               soil_class    = NULL,
                               lookup        = NULL,
                               cpus          = 1) {
  if (!dir.exists('class_maps/')) {
    dir.create('class_maps/', showWarnings = F)
  }

  class_dir <- file.path(getwd(), 'class_maps')

  suppressWarnings(
    probs_list <- if(is.null(soil_class)) {
      raster::unstack(tallied_probs)
      } else {
        raster::unstack(raster::subset(tallied_probs, soil_class))
        })

  if(is.null(soil_class)) {
  lapply(seq_along(probs_list), function(i) {
    names(probs_list[[i]]) <- as.character(lookup[i, 2])
  })
  } else {
     names(probs_list) <- soil_class
    }

  prob_tifs <- lapply(probs_list, function(x) {
    writeRaster(x,
                filename  = file.path(class_dir, paste0(names(x), "_probability.tif")),
                format    = "GTiff",
                NAflag    = -9999,
                datatype  = 'FLT4S',
                overwrite = TRUE)
  })

  assign('class_maps', prob_tifs, envir = .GlobalEnv)
}
