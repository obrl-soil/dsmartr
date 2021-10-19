#' Prepare dsmartr polygons
#'
#' Prepares soil map inputs for use in
#' \code{\link[dsmartr:iterate]{dsmartr::iterate()}}.
#' @param src_map An sfc_POLYGON/MULTIPOLYGON or SpatialPolygonsDataFrame object
#'   representing the soil map to be disaggregated. Format requirements:
#' \itemize{
#' \item{one row per polygon (data is wide-formatted)}
#' \item{A numeric unique ID field for polygons}
#' \item{1-n character columns for soil classes named CLASS_1 to CLASS_n}
#' \item{1-n numeric columns for soil class percentages named PERC_1 to PERC_n.
#' PERC_1 must relate to CLASS_1, etc.}}
#'
#' Other columns may exist in the object; they will be ignored.
#' @param covariates RasterStack or RasterBrick; environmental covariate
#'   datasets.
#' @param id_field String; name of unique identifier field in \code{src_map}.
#' @param sample_method String; choice of flat rate per polygon or
#'   area-proportional rate.
#' @param sample_rate Integer; Number of samples per polygon.
#' @param rate_floor Integer; desired minimum number of samples per polygon.
#'   Optional; used with \code{sample_method = 'area_p'}. Defaults to 2x the
#'   number of soil classes on a polygon, switch off with samp_floor = 0 (not
#'   recommended).
#' @param rate_ceiling Integer; desired maximum number of samples per polygon.
#'   Optional; use with \code{sample_method = 'area_p'}. Only applies to
#'   polygons with a single class component and (effectively) a large area.
#' @return A data frame holding polygon input attributes and four new attribute
#'   columns:
#' \itemize{
#'   \item{\code{area_sqkm}: Polygon area in square kilometers, by
#'   \code{\link[sf:geos_measures]{sf::st_area()}}.}
#'   \item{\code{n_soils}: The number of soil classes within the map unit.}
#'   \item{\code{n_samples}: The number of environmental covariate point samples
#'   that will be taken on each model iteration.}
#'   \item{\code{intersecting_cells}: Raster cell index numbers for any cell
#'   whose center falls within the polygon boundary.}
#'   }
#' Outputs are also written to disk.
#' @note \itemize{
#'   \item{The output of this function is a required input for
#'   \code{\link[dsmartr:iterate]{dsmartr::iterate()}}}.
#'   \item{Covariate data should be in a projected CRS with defined units. While
#'   the function will run with lat/long data, polygon area calculations may be
#'   dangerously inaccurate. Vector inputs will be transformed to match the
#'   covariate CRS.}
#'   \item{The \code{intersecting_cells} attribute field is a list-column, so
#'   the returned object cannot be written to e.g. csv format.}
#' \item{This function runs faster with a RasterBrick than a Stack.}}
#' @examples \dontrun{
#' data('heronvale_soilmap')
#' data('heronvale_covariates')
#'
#' # flat rate
#' pr_flat <- prep_polygons(src_map = heronvale_soilmap,
#'                          covariates = heronvale_covariates,
#'                          id_field = 'POLY_NO', sample_method = 'flat',
#'                          sample_rate = 6)
#'
#' # area_proportional rate with floor
#' pr_ap <- prep_polygons(src_map = heronvale_soilmap,
#'                        covariates = heronvale_covariates,
#'                        id_field = 'POLY_NO', sample_method = 'area_p',
#'                        sample_rate = 20, rate_floor = 6)
#' }
#' @importFrom methods as
#' @importFrom sf st_area st_set_geometry st_transform
#' @importFrom raster crs
#' @importFrom units as_units
#' @export
#'
prep_polygons <- function(src_map       = NULL,
                          covariates    = NULL,
                          id_field      = NULL,
                          sample_method = c('flat', 'area_p'),
                          sample_rate   = NULL,
                          rate_floor    = NULL,
                          rate_ceiling  = NULL) {

  sample_method <- match.arg(sample_method)

  if (!dir.exists(file.path(getwd(), 'inputs'))) {
    dir.create(file.path(getwd(), 'inputs'), showWarnings = FALSE)
  }
  in_dir <- file.path(getwd(), 'inputs')

  # coerce sp polygons to sf
  if(!inherits(src_map, 'sf')) {
    src_map <- sf::st_as_sf(src_map)
  }

  fct <- sapply(src_map, is.factor)
  src_map[fct] <- lapply(src_map[fct], as.character)
  src_prepped <-
    sf::st_transform(src_map, crs = raster::crs(covariates, asText = TRUE))

  # polygon area in sq km
  src_prepped$area_sqkm <- sf::st_area(src_prepped)
  units(src_prepped$area_sqkm) <- units::as_units('km^2')

  # number of soil classes on polygon
  src_prepped$n_soils <- apply(sf::st_set_geometry(src_prepped, NULL),
                               MARGIN = 1, FUN = function(row) {
                                 sum(!is.na(row) &
                                       grepl('CLASS', names(row)) == TRUE)
                                 })
  # intended sample number
  if(sample_method == 'flat') {
    src_prepped$n_samples <- sample_rate
  }

  if(sample_method == 'area_p') {
    src_prepped$temparea <- as.numeric(src_prepped$area_sqkm) # apply coercion :|
    src_prepped$n_samples <-
      apply(sf::st_set_geometry(src_prepped, NULL), MARGIN = 1, function(row) {
        p_nsoils   <- as.integer(row[['n_soils']])
        p_area     <- as.numeric(row[['temparea']])
        samp_area  <- ceiling(p_area * sample_rate)
        samp_floor <- if(is.null(rate_floor)) { p_nsoils * 2 } else { rate_floor }
        if(p_nsoils == 1L & !is.null(rate_ceiling)) {
          as.integer(max(min(rate_ceiling, samp_area), samp_floor))
        } else {
          as.integer(max(samp_area, samp_floor))
        }
      })
    src_prepped$temparea <- NULL
    }

  # list cell numbers that intersect the polygon
  # (these will be subset randomly on each model run)
  src_prepped$intersecting_cells <- strict_cfp(src_map    = src_prepped,
                                               covariates = covariates)

  ## old way - buggy, some problem in raster pkg, src unclear at Jan 17
  #intersecting_cells <- map(src_split, function(g) {
  #             poly_sp <- as(g, 'Spatial')
  #             cells   <- as.integer(unlist(cellFromPolygon(covariates, poly_sp),
  #                                          use.names = FALSE))
  #             cells   <- if (length(cells) == 0) { NA } else { cells }
  #           })

  # can't index a list-column
  src_prepped <- src_prepped[which(!is.na(src_prepped$intersecting_cells)), ]
  src_prepped <- sf::st_set_geometry(src_prepped, NULL)

  IDstr_r <- if (nrow(src_map) > nrow(src_prepped)) {
    removed   <- setdiff(src_map[[id_field]], src_prepped[[id_field]])
    n_removed <- length(removed)
    IDstr     <- toString(unlist(removed))
    message(paste0('Polygons with ID(s) ', IDstr, ' excluded - no cell intersections.'))
    IDstr
  } else {
    NA
  }

  saveRDS(src_prepped, file.path(in_dir, 'prepped_map.rds'))
  src_prepped

}

#' Prepare dsmartr points
#'
#' Prepares dsmartr point inputs for use in
#' \code{\link[dsmartr:iterate]{dsmartr::iterate()}}.
#' @param known_points sfc_POINT, SpatialPointsDataFrame or Data Frame object;
#'   represents locations where the soil class has been directly observed.
#'   Supply with spatial location (spatial data or x and y attributes), unique
#'   numeric ID, and soil class.
#' @param soil_id String; name of column in \code{known_points} holding soil
#'   class data.
#' @param x_coords String; name of column in \code{known_points} holding
#'   x-coordinate data. Only needed for non-spatial inputs.
#' @param y_coords String; name of column holding y-coordinate data. Only needed
#'   for non-spatial inputs.
#' @param covariates RasterStack or RasterBrick; environmental covariate data.
#' @return A data frame with two attributes: soil class code and corresponding
#'   raster cell index number.
#' @note The output of this function is an optional input for
#'   \code{\link{iterate}}. If \code{known_points} is supplied as a data frame,
#'   the x and y coordinates must match the covariate crs.
#' @examples \dontrun{
#' load('heronvale_known_sites')
#' load('heronvale_soilmap')
#' load('heronvale_covariates')
#'
#' # data frame input:
#' pr_pts <- prep_points(known_points = heronvale_known_sites, soil_id = 'CLASS',
#' x_coords = 'x', y_coords = 'y', covariates = heronvale_covariates)
#' }
#' @importFrom methods as
#' @importFrom raster cellFromXY crs extract
#' @importFrom sp spTransform
#' @importFrom sf st_as_sf st_set_geometry st_transform
#' @export
#'
prep_points <- function(known_points = NULL, soil_id  = NULL,
                        x_coords     = NULL, y_coords = NULL,
                        covariates   = NULL) {

  if (!dir.exists(file.path(getwd(), 'inputs'))) {
    dir.create(file.path(getwd(), 'inputs'), showWarnings = FALSE)
  }
  in_dir <- file.path(getwd(), 'inputs')

  fct <- sapply(known_points, is.factor)
  known_points[fct] <- lapply(known_points[fct], as.character)

  kpp <- if(inherits(known_points, 'sf')) {
    known_points <-
      sf::st_transform(known_points, crs = raster::crs(covariates,
                                                       asText = TRUE))
    cellnos <- raster::cellFromXY(covariates, xy = as(known_points, 'Spatial'))
    known_points$CELL <- cellnos
    known_points  <- sf::st_set_geometry(known_points, NULL)
    names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
    known_points <- known_points[ , c('CLASS', 'CELL')]

    } else if (inherits(known_points, 'Spatial')) {
       known_points <-
         sp::spTransform(known_points, CRSobj = raster::crs(covariates))
       cellnos <- raster::cellFromXY(covariates, xy = known_points)
       known_points <- known_points@data
       names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
       known_points$CELL  <- cellnos
       known_points <- known_points[ , c('CLASS', 'CELL')]

       } else if(inherits(known_points, 'data.frame')) {
         known_points <-
           sf::st_as_sf(known_points, coords = c(x_coords, y_coords),
                        crs = raster::crs(covariates, asText = TRUE))
         cellnos <- raster::cellFromXY(covariates,
                                       xy = as(known_points, 'Spatial'))
         names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
         known_points$CELL <- cellnos
         known_points <- known_points[ , c('CLASS', 'CELL')]

       } else { stop('Please supply a valid input file') }

  saveRDS(kpp, file.path(in_dir, 'prepped_points.rds'))
  kpp
}
