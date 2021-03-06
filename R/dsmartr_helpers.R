#' get items
#'
#' Returns cell values from selected columns in an attribute table as an atomic
#' vector, with NA values excluded.
#' @keywords internal
#' @param input One row of a data frame, as vector
#' @param selector String; a column name or name stub.
#' @return An atomic vector of non-NA values in the selected columns, by order
#'   of occurrence.
#' @note This is a helper function and not widely applicable. It is expected to
#'   be used by-row on wide-formatted spatial attribute data.
#' @examples \dontrun{
#' load('heronvale_soilmap')
#' sm <- sf::st_set_geometry(heronvale_soilmap, NULL)
#' percs_1 <- n_things(sm[1, ], 'PERC')
#' }
#' @importFrom stats na.omit
#'
n_things <- function(input = NULL, selector = NULL) {
  if(length(grep(selector, names(input))) == 0) {
    stop("Error: Selector not found within input's column names")
  } else {
  output <- as.vector(na.omit(unlist(input[c(grep(selector, names(input)))])))
  if(length(output) == 0) { NA } else { output }}
}

#' in range
#'
#' Returns true if number is within range
#' @param value Numeric; number to check
#' @param lower Numeric; lower limit of range
#' @param upper Numeric; upper limit of range
#' @param strict Logical; < or <= ?
#' @return TRUE if value is within range
#' @examples
#' in_range(5, 3, 6)
#' in_range(5, 3, 5, FALSE)
#' in_range(5, 3, 5, TRUE)
#' @export
#'
in_range <- function(value = NULL, lower = NULL, upper = NULL, strict = TRUE) {
  stopifnot(lower < upper)
  if(strict) {
    (value - lower) * (upper - value) > 0
  } else {
    (value - lower) * (upper - value) >= 0
  }
}

#' Check dsmartr input map attributes
#'
#' Checks dsmartr input map attributes for potential data quality issues.
#' @param src_map An sfc_POLYGON/MULTIPOLYGON or SpatialPolygonsDataFrame
#'   object, see \code{\link[dsmartr:prep_polygons]{dsmartr::prep_polygons()}}.
#' @param id_field String; name of unique feature identifier field.
#' @param cs String; stub of attribute column names holding soil class data
#' @param ps String; stub of attribute column names folding soil percentage data
#' @details This function highlights some issues that may cause dsmartr to
#'   either fail or behave unpredictably. Usually these are side-effects of
#'   prepration work carried out in other software. Sometimes they are the
#'   result of errors in the source data. Note that only attributes are checked;
#'   geometry checking techniques are already available in sf and other spatial
#'   packages. At present, four problems are checked: whether data is missing
#'   (e.g. a polygon class without a matching percentage), whether any
#'   percentage values are 0\%, whether total percentage (by polygon) is not
#'   100\%, and whether any IDs are duplicated.
#' @return An sf object with appended logical columns indicating presence of
#'   possible data faults.
#' @examples \dontrun{
#' load("heronvale_soilmap")
#' checked_map <- check_attributes(src_map  = heronvale_soilmap,
#'                                 id_field = 'POLY_NO',
#'                                 cs = 'CLASS',
#'                                 ps = 'PERC')
#'
#' # to write to file without data loss (some formats
#' # can't store boolean), do this beforehand:
#' checked_map <- dplyr::mutate_if(checked_map, is.logical, as.character)
#' }
#' @importFrom sf st_as_sf st_set_geometry
#' @export
#'
check_attributes <- function(src_map  = NULL, id_field = NULL,
                             cs = NULL, ps = NULL) {

  # coerce to sf
  if(is(src_map, 'sf') == FALSE) {
    src_map <- sf::st_as_sf(src_map)
  }

  sm <- sf::st_set_geometry(src_map, NULL)

  ### put problems in new columns for easy ID

  src_map$missing_data <- apply(sm, MARGIN = 1, FUN = function(row) {
    n_classes <- length(n_things(row, cs))
    n_percs   <- length(n_things(row, ps))
    if(n_classes != n_percs) { TRUE } else { FALSE }
    })

  src_map$zero_percs <- apply(sm, MARGIN = 1, FUN = function(row) {
    poly_percs <- as.numeric(n_things(row, ps))
    if(any(poly_percs == 0)) { TRUE } else { FALSE }
    })

  src_map$problem_percs <- apply(sm, MARGIN = 1, FUN = function(row) {
               total <- sum(as.numeric(n_things(row, ps)))
               if(total != 100) { TRUE } else { FALSE }
               })

  src_map$duplicate_ids <- duplicated(sm[[id_field]])

  sf::st_as_sf(src_map) # re-orders cols so geom at end
}

#' generate masks from samples
#'
#' Generates a rasterlayer for each dsmart iteration highlighting areas where
#' one or more covariates has a value that is out of range of the sample used to
#' model soil distribution.
#' @param samples List; POINT_sfc objects created by
#'   \code{\link[dsmartr:iterate]{dsmartr::iterate()}} when option write_samples
#'   = TRUE
#' @param covariates RasterStack or RasterBrick; environmental covariate
#'   datasets.
#' @param tolerance Integer; the number of out-of-range covariates allowed.
#'   Defaults to 0, that is, if any covariate is out of range on a given pixel,
#'   that pixel will be masked.
#' @param cpus Integer; number of processors to use in parallel.
#' @details For each model iteration this function generates a raster that can
#'   be used to mask dsmartr predictions where covariate values are out of range
#'   of the sample of cells used in that iteration. The masks can also be
#'   combined to produce an overall mask, applicable to final products. This is
#'   experimental; the idea is to avoid making predictions in locations where
#'   the model won't perform well.
#' @return a list of rasters
#' @examples \dontrun{
#' # run dsmartr::iterate() with the example data, then:
#' sample_list <- lapply(list.files(file.path(getwd(), 'iterations', 'samples'),
#' pattern = '\\.gpkg', full.names = TRUE), read_sf)
#' masks <- prediction_masks(samples = sample_list,
#'                           covariates = heronvale_covariates,
#'                           cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # less strict: up to 3 out of range covariates are allowed
#' masks <- prediction_masks(samples = sample_list,
#'                           covariates = heronvale_covariates,
#'                           tolerance = 3,
#'                           cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # a 'final product' mask, where outputs are masked only if the mode of that
#' # cell is 'do not predict'.
#' all_masks <- stack(masks)
#' modal_mask <- calc(all_masks, function(cell) ifelse(raster::modal(cell) == 1,
#' NA, 0))
#' # the above could be applied to the most-likely soils map after running
#' # dsmartr::most_likely:
#' masked_m1 <- most_likely_soil[[1]][[1]] + modal_mask
#' }
#' @importFrom raster calc writeRaster
#' @importFrom sf st_set_geometry
#' @export
#'
prediction_masks <- function(samples = NULL, covariates = NULL,
                             tolerance = 0L, cpus = 1) {

  if (!dir.exists(file.path(getwd(), 'iterations', 'masks'))) {
    dir.create(file.path(getwd(), 'iterations', 'masks'),
               showWarnings = FALSE, recursive = TRUE)
  }
  drmsk <- file.path(getwd(), 'iterations', 'masks')
  if(length(samples) == 0) {
    stop('Samples missing.')
  }
  message(paste0(Sys.time(), ': dsmartr prediction mask creation in progress...'))
  pb <- txtProgressBar(min = 0, max = length(samples), style = 3)

  all_masks <- mapply(FUN = function(m, n) {
    sx <- st_set_geometry(m, NULL)
    rangex <- apply(sx[, c(3:ncol(sx))], MARGIN = 2,
                    FUN = function(x) { range(x, na.rm = TRUE) } )
    sx_mins <- rangex[1,]
    sx_maxs <- rangex[2,]

    beginCluster(n = cpus)
    assign('sx_mins', sx_mins, envir = parent.frame())
    assign('sx_maxs', sx_maxs, envir = parent.frame())
    pred_mask <-
      clusterR(x    = covariates,
               fun  = calc,
               args = list(fun = function(cell) {
                 in_or_out <- in_range(cell, sx_mins, sx_maxs, strict = FALSE)
                 # handle no data areas
                 if(all(is.na(in_or_out))) {
                   return(NA_integer_)
                 }
                 if(sum(in_or_out == FALSE, na.rm = TRUE) > tolerance) {
                     1L
                   } else {
                     0L
                     }
                 }),
               filename  = file.path(drmsk, paste0('pred_mask_', n, '.tif')),
               NAflag    = -9999,
               datatype  = 'INT2S',
               overwrite = TRUE)
    setTxtProgressBar(pb, n)
    endCluster()
    pred_mask
    },
    m = samples,
    n = seq_along(samples))

  close(pb)

  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 drmsk))
  all_masks
}

#' get a list of cell numbers intersecting a polygon
#'
#' Determines which covariate cells intersect with a given map polygon.
#' @param src_map An sfc_POLYGON/MULTIPOLYGON or SpatialPolygonsDataFrame object
#'   representing the soil map to be disaggregated. Attributes are not required.
#' @param covariates RasterStack or RasterBrick; environmental covariate
#'   datasets.
#' @details Given a polygon layer and spatially intersecting raster object,
#'   determines which raster cells intersect which polygons. This function is
#'   primarily intended for use in \code{\link{prep_polygons}} but can be used
#'   as a standalone.
#' @return a list of atomic vectors containing raster cell numbers.
#' @examples \dontrun{
#' data('heronvale_soilmap')
#' data('heronvale_covariates')
#'
#' # just get the list object
#' intersecting_cells <- strict_cfp(src_map = heronvale_soilmap,
#'                                  covariates = heronvale_covariates)
#'
#' # to check that this output is complete and non overlapping
#' r_blank <- raster(heronvale_covariates)
#' r_blank[] <- 1:ncell(r_blank)
#' dup_chk <- table(unlist(intersecting_cells)) # all shld be freq = 1
#' values(r_blank)[!(values(r_blank) %in% as.data.frame(dup_chk)$Var1)] <- NA
#' plot(r_blank) # only cells underlying the polygons are kept
#'
#' # or place it in a new list-column for later use
#' heronvale_soilmap$int_cells <- strict_cfp(src_map = heronvale_soilmap,
#'                                           covariates = heronvale_covariates)
#'
#' # as above, tidyverse style
#' library(dplyr)
#' heronvale_soilmap <- heronvale_soilmap %>%
#'   mutate(int_cells = strict_cfp(., heronvale_covariates))
#'
#' }
#' @importFrom raster alignExtent crop extent intersect ncell raster
#'   rasterToPoints
#' @importFrom sf st_geometry st_intersection st_as_sf
#' @export
#'
strict_cfp <- function(src_map = NULL, covariates = NULL) {

  # init empty raster
  r_blank <- raster::raster(covariates)
  # populate with cell index values
  r_blank[] <- seq.int(ncell(r_blank))

  inps <- split(src_map, seq(nrow(src_map)))

  # lapply and map produce identical output in identical times here
  lapply(inps, FUN = function(p) {
    p_extent <- alignExtent(extent(p), r_blank, snap = 'near')
    # catch polygons off the edge of a raster
    r_blank_p <- try(raster::crop(r_blank, p_extent)) # src of rgeos dep
    if(inherits(r_blank_p, 'try-error')) {
      return(NA_integer_)
    }
    blank_pts <- raster::rasterToPoints(r_blank_p, spatial = TRUE)
    blank_pts <- sf::st_as_sf(blank_pts)
    out <- suppressWarnings(sf::st_intersection(blank_pts, sf::st_geometry(p)))
    as.integer(out$layer)
   })
}
