#' Extract most-likely-soil maps
#'
#' This function takes the outputs of \code{\link{collate}} and extracts the top
#' \code{n} most-likely soil maps. Optionally, their respective probability
#' surfaces are also extracted.
#' @param dsmart_preds RasterBrick; output \code{dsmartr_predictions} from
#'   \code{\link[dsmartr:collate]{dsmartr::collate()}}.
#' @param dsmart_probs RasterBrick; output \code{dsmartr_probabilities} from
#'   \code{\link[dsmartr:collate]{dsmartr::collate()}}; optional.
#' @param n_maps integer; the number of most-probable maps desired.
#' @return A list containing \code{n_maps} most-likely-soils maps, from most to
#'   least probable. Optionally; \code{n_maps} probability surfaces associated
#'   with \code{most_likely_maps} follow.
#' All outputs are written to disk as GeoTIFFs.
#' @examples \dontrun{
#' # run collate() example code, then:
#' most_likely_soil <-
#'  most_likely(dsmart_preds = collated[['dsmartr_predictions']],
#'              dsmart_probs = collated[['dsmartr_probabilities']],
#'              n_maps = 2)
#' }
#' @importFrom raster unstack ratify writeRaster
#' @export
#'
most_likely <- function(dsmart_preds = NULL,
                        dsmart_probs = NULL,
                        n_maps       = 2L) {
  if (!dir.exists(file.path(getwd(), 'most_likely_maps'))) {
    dir.create(file.path(getwd(), 'most_likely_maps'), showWarnings = FALSE)
  }
  mp_dir <- file.path(getwd(), 'most_likely_maps')

  message(paste0(Sys.time(), ': dsmartr prediction map unstacking in progress...'))

  if (!is.null(dsmart_preds)) {
    suppressWarnings(map_list <- raster::unstack(dsmart_preds[[1:n_maps]]))
    names(map_list) <- paste0('most_likely_', 1:n_maps)
  } else {
    stop('dsmart_preds raster stack has not been specified')
  }

  pb <- txtProgressBar(min = 0, max = n_maps, style = 3)
  most_likely_maps <- mapply(FUN = function(x, i) {
    ml <- writeRaster(x,
                filename  = file.path(mp_dir, paste0('most_likely_', i ,'.tif')),
                format    = 'GTiff',
                NAflag    = -9999,
                datatype  = 'INT2S',
                overwrite = TRUE)
    setTxtProgressBar(pb, i)
    ml
  },
  x = map_list, i = seq_along(1:n_maps)
  )
  close(pb)

  ### optional probability surfaces
  if(!is.null(dsmart_probs)) {
    message(paste0(Sys.time(), ': dsmartr probability surface unstacking in progress...'))
    suppressWarnings(probmap_list <- raster::unstack(dsmart_probs[[1:n_maps]]))
    names(probmap_list) <- paste0('most_likely_prob_', 1:n_maps)

    pb <- txtProgressBar(min = 0, max = n_maps, style = 3)
    most_likely_ps <- mapply(FUN = function(x, i) {
      mp <- writeRaster(x,
                  filename  = file.path(mp_dir, paste0('most_likely_prob_', i, '.tif')),
                  format    = 'GTiff',
                  datatype  = 'FLT4S',
                  NAflag    = -9999,
                  overwrite = TRUE)
      setTxtProgressBar(pb, i)
      mp
    },
    x = probmap_list,
    i = seq_along(1:n_maps)
    )
    close(pb)
    message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                   file.path(getwd(), 'most_likely_maps')))
  } else {
    message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                   file.path(getwd(), 'most_likely_maps')))
  }

  unstacked_maps <- if(is.null(dsmart_probs)) {
    most_likely_maps
  } else {
    c(most_likely_maps, most_likely_ps)
  }

}

#' Produce soil class probability surfaces
#'
#' Generates a soil class probability map for any or all input soil classes.
#' Requires output of
#' \code{\link{collate}}.
#' @param tallied_probs RasterBrick; \code{tallied_probabilities} output by
#'   \code{\link{collate}}.
#' @param soil_class String; Soil class(es) of interest, if only particular maps
#'   are desired.
#' @return \code{class_maps}: List of RasterLayers; probability surfaces for
#'   each soil class. All outputs are written to disk as GeoTIFF.
#' @examples \dontrun{
#' # run collate() example code, then
#'
#' # all classes:
#' classmaps <- class_maps(tallied_probs = collated[['tallied_probabilities']])
#'
#' # just map two classes of interest:
#' two_classmaps <- class_maps(tallied_probs = collated[['tallied_probabilities']],
#' soil_class = c('An', 'Wr'))}
#' @importFrom raster unstack writeRaster
#' @export
class_maps <- function(tallied_probs = NULL,
                       soil_class    = NULL) {
  if (!dir.exists(file.path(getwd(), 'class_maps'))) {
    dir.create(file.path(getwd(), 'class_maps'), showWarnings = FALSE)
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
    names(probs_list) <- names(tallied_probs)
    } else {
      names(probs_list) <- soil_class
    }

  pb <- txtProgressBar(min = 0, max = length(probs_list), style = 3)
  class_ps <- lapply(seq_along(probs_list), function(x) {
    cm <- writeRaster(probs_list[[x]],
                      filename  = file.path(class_dir,
                                            paste0(names(probs_list)[[x]], "_probability.tif")),
                      format    = "GTiff",
                      NAflag    = -9999,
                      datatype  = 'FLT4S',
                      overwrite = TRUE)
    setTxtProgressBar(pb, x)
    cm
  })
  close(pb)
  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'class_maps')))
  class_ps
}
