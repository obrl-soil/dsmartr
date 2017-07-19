#' Extract most-likely-soil maps
#'
#' This function takes the outputs of \code{\link{dsmartr_collate}} and extracts the top \code{n}
#' most-likely soil maps. Optionally, their respective probability surfaces are also extracted.
#' @param dsmart_preds RasterBrick; output \code{dsmart_predictions} from
#'   \code{\link{dsmartr_collate}}.
#' @param dsmart_probs RasterBrick; output \code{dsmart_probabilities} from
#'   \code{\link{dsmartr_collate}}; optional.
#' @param n_maps integer; the number of most-probable maps desired.
#' @return A list containing:
#' \itemize{
#'   \item{\code{most_likely_maps}: A list of \code{n_maps} most-likely-soils maps, from most to least
#'   probable.}
#'   \item{\code{most_likely_ps}: Optional; A list of \code{n_maps} probability surfaces
#'   associated with \code{most_likely_maps}.}}
#' All outputs are written to disk as GeoTIFFs.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data and then:
#' most_likely <- dsmartr_most_likely(dsmart_preds = collated[['dsmartr_predictions']],
#' dsmart_probs = collated[['dsmartr_probabilities']], n_maps = 2)}
#' @importFrom raster unstack ratify writeRaster
#' @export
dsmartr_most_likely <- function(dsmart_preds = NULL,
                                  dsmart_probs = NULL,
                                  n_maps       = 2) {
  if (!dir.exists('most_probable_maps/')) {
    dir.create('most_probable_maps/', showWarnings = F)
  }
  mp_dir <- file.path(getwd(), 'most_probable_maps')

  message(paste0(Sys.time(), ': dsmartr prediction map unstacking in progress...'))

  if (!is.null(dsmart_preds)) {
    suppressWarnings(map_list <- raster::unstack(dsmart_preds[[1:n_maps]]))
  } else {
    stop('dsmart_preds raster stack has not been specified')
  }

  lapply(map_list, function(i) ratify(i))

  most_likely_maps <- mapply(FUN = function(x, i) {
    writeRaster(x,
                filename  = file.path(mp_dir, paste0('mostlikely_', i ,'.tif')),
                format    = 'GTiff',
                NAflag    = -9999,
                datatype  = 'INT2S',
                overwrite = TRUE)
  },
  x = map_list, i = seq_along(1:n_maps)
  )

  ### optional probability surfaces
  if(!is.null(dsmart_probs)) {
    message(paste0(Sys.time(), ': dsmartr probability surface unstacking in progress...'))
    suppressWarnings(probmap_list <- raster::unstack(dsmart_probs[[1:n_maps]]))

    most_likely_ps <- mapply(FUN = function(x, i) {
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
    message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                   file.path(getwd(), 'most_probable_maps')))
  } else {
    message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                   file.path(getwd(), 'most_probable_maps')))
  }

  unstacked_maps <- if(is.null(dsmart_probs)) {
    most_likely_maps
  } else {
    list("most_likely_maps" = most_likely_maps,
         "most_likely_ps"   = most_likely_ps ) # ps 'probability surface'
  }

}

#' Produce soil class probability surfaces
#'
#' Generates a soil class probability map for any or all input soil classes. Requires output of
#' \code{\link{dsmartr_collate}}.
#' @param tallied_probs RasterBrick; \code{tallied_probabilities} output by
#'   \code{\link{dsmartr_collate}}.
#' @param soil_class String; Soil class(es) of interest, if only particular maps are desired.
#' @param lookup Data Frame; contains raster values and corresponding soil class labels. See
#'   \code{\link{dsmartr_collate}}.
#' @param cpus Integer; number of processors to use in parallel.
#' @return \code{class_maps}: List of RasterLayers; probability surfaces for each soil class. All
#'   outputs are written to disk as GeoTIFF.
#' @examples \dontrun{
#' # run dsmartr_collate() with the example data and then
#'
#' # all classes:
#' classmaps <- dsmartr_class_maps(tallied_probs = collated[['tallied_probabilities']],
#'  lookup = LUT)
#'
#' # just map two classes of interest:
#' classmaps <- dsmartr_class_maps(tallied_probs = collated[['tallied_probabilities']],
#' soil_class = c('An', 'Wr'), lookup = LUT)}
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

  message(paste0(Sys.time(), ': dsmartr class map production in progress...'))
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

  class_ps <- lapply(probs_list, function(x) {
    writeRaster(x,
                filename  = file.path(class_dir, paste0(names(x), "_probability.tif")),
                format    = "GTiff",
                NAflag    = -9999,
                datatype  = 'FLT4S',
                overwrite = TRUE)
  })

  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'class_maps')))
  class_ps
}
