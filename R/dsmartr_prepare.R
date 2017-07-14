#' Prepare dsmartr polygons
#'
#' Prepares dsmartr polygon inputs for use in \code{\link{dsmartr_iterate}}.
#' @param src_map An sfc_POLYGON, sfc_MULTIPOLYGON, or SpatialPolygonsDataFrame object representing
#' the soil map to be disaggregated. Format requirements:
#' \itemize{
#' \item{one row per polygon (data is wide-formatted)}
#' \item{A numeric unique ID field for polygons}
#' \item{1-n character columns for soil classes named CLASS_1 to CLASS_n}
#' \item{1-n numeric columns for soil class percentages named PERC_1 to PERC_n. PERC_1 must relate
#' to CLASS_1, etc.}
#' }
#' Other columns may exist in the object; they will be ignored.
#' @param covariates RasterStack or RasterBrick; environmental covariate data.
#' @param id_field String; name of unique identifier field in \code{src_map}.
#' @param sample_method String; choice of flat rate per polygon or area-proportional rate.
#' @param flat_rate Integer; Number of samples per polygon; use with \code{sample_method = 'flat'}.
#' @param area_rate Integer; desired number of samples per square kilometre; use with
#' \code{sample_method = 'area_p'}.
#' @param floor Integer; desired minimum number of samples per polygon. Optional; use with
#' \code{sample_method = 'area_p'}. Defaults to the number of soil classes on a polygon.
#' @return A data frame holding polygon input attributes and four new attribute columns:
#' \itemize{
#'   \item{\code{area_sqkm}: Polygon area in square kilometers, by \code{\link[sf]{st_area}}.}
#'   \item{\code{n_classes}: The number of soil classes within the map unit.}
#'   \item{\code{n_samples}: The number of environmental covariate point samples that will be taken
#'   on each model iteration.}
#'   \item{\code{intersecting_cells}: Raster cell index numbers for any cell whose center falls
#'   within the polygon boundary.}
#'   }
#' Outputs are also written to disk.
#' @note \itemize{
#' \item{The output of this function is a required input for \code{\link{dsmartr_iterate}}}.
#' \item{Covariate data must be in a projected CRS with defined units (not lat/long). Vector inputs
#' will be transformed to match the covariate CRS.}
#' \item{The \code{intersecting_cells} attribute field is a list-column, so the returned  object
#' cannot be written to e.g. csv format.}
#' \item{This function runs faster with a RasterBrick than a Stack.}}
#' @examples \dontrun{
#' load("heronvale_soilmap")
#' load("heronvale_covariates")
#'
#' # flat rate
#' dsmartr_prep_polygons(src_map = heronvale_soilmap, covariates = heronvale_covariates,
#' id_field = 'POLY_NO', sample_method = 'flat', flat_rate = 20)
#'
#' # area_proportional rate with floor
#' dsmartr_prep_polygons(src_map = heronvale_soilmap, covariates = heronvale_covariates,
#' id_field = 'POLY_NO', sample_method = 'area_p', area_rate = 20, floor = 10)}
#' @importFrom dplyr filter mutate mutate_if
#' @importFrom purrr map map_int
#' @importFrom stats setNames
#' @importFrom raster crs cellFromPolygon
#' @export
dsmartr_prep_polygons <- function(src_map       = NULL,
                                  covariates    = NULL,
                                  id_field      = NULL,
                                  sample_method = c('flat', 'area_p'),
                                  flat_rate     = NULL,
                                  area_rate     = NULL,
                                  floor         = NULL) {

  if (!dir.exists('inputs')) {
    dir.create('inputs', showWarnings = F)
  }
  in_dir <- file.path(getwd(), 'inputs')

  # coerce polygons to sf
  src_map <- if(is(src_map, 'sf') == FALSE) {
    st_as_sf(src_map)
  } else {
    src_map
  }

  src_map <- mutate_if(src_map, is.factor, as.character)
  src_map <- st_transform(src_map, crs = crs(covariates, asText = TRUE))

  # polygon area in sq km
  src_prepped <- mutate(src_map, area_sqkm = st_area(src_map))
  units(src_prepped$area_sqkm) <- with(units::ud_units, km^2)
  src_split <- split(src_prepped, 1:nrow(src_prepped))

  # number of soil classes on polygon
  n_soils <- map_int(src_split, function(ncl) {
               sum(!is.na(ncl) & grepl('CLASS', names(ncl)) == T)
             })

  src_prepped$n_soils <- as.vector(n_soils)
  src_split <- split(src_prepped, 1:nrow(src_prepped))

  # intended sample number
  n_samples <- if(sample_method == 'flat') {
    flat_rate
    } else if (sample_method == 'area_p') {
      map_int(src_split, function(ns) {
               p_percs   <- as.numeric(n_things(ns, 'PERC'))
               p_nsoils  <- as.integer(as.data.frame(ns)[, 'n_soils'])
               p_area    <- as.numeric(as.data.frame(ns)[, 'area_sqkm'])
               gap       <- max(p_percs) / min(p_percs)
               samp_min  <- ceiling(max(1, gap * p_nsoils))
               samp_area <- ceiling(p_area * area_rate)
               floor     <- if(is.null(floor)) { p_nsoils } else { floor }
               samp_n    <- as.integer(max(samp_min, samp_area, floor))
             } )
    } else {
      stop('Please provide a valid sample_method parameter. Options are \'flat\', \'area_p\'.')
    }

  src_prepped$n_samples <- as.vector(n_samples)
  src_split <- split(src_prepped, 1:nrow(src_prepped))

  # list cell numbers that intersect the polygon
  #(these will be subset randomly on each model run)
  intersecting_cells <- map(src_split, function(g) {
               poly_sp <- as(g, 'Spatial')
               cells   <- as.integer(unlist(cellFromPolygon(na.omit(covariates), poly_sp)))
               cells   <- if (length(cells) == 0) { NA } else { cells }
             })
  intersecting_cells <- setNames(intersecting_cells, NULL)

  src_prepped$intersecting_cells <- intersecting_cells #list-column

  src_prepped <- filter(src_prepped, !(is.na(intersecting_cells))) # can't index a list-column

  src_prepped <- st_set_geometry(src_prepped, NULL)

  IDstr_r <- if (nrow(src_map) > nrow(src_prepped)) {
    removed   <- setdiff(as.data.frame(src_map)[ , id_field],
                         src_prepped[ , id_field])
    n_removed <- nrow(removed)
    IDstr     <- toString(unlist(removed))
    message(paste0('Polygons with ID(s) ', IDstr, ' excluded - no cell intersections.'))
    IDstr
  } else {
    NA
  }

  saveRDS(src_prepped, file.path(in_dir, 'prepped_map.rds'))
  return(src_prepped)

}

#' Prepare dsmartr points
#'
#' Prepares dsmartr point inputs for use in \code{\link{dsmartr_iterate}}.
#' @param known_points sfc_POINT, SpatialPointsDataFrame or Data Frame object; represents locations
#'  where the soil class has been directly observed. Supply with spatial location (spatial data or
#'  x and y attributes), unique numeric ID, and soil class.
#' @param soil_id String; name of column in \code{known_points} holding soil class data.
#' @param x_coords String; name of column in \code{known_points} holding x-coordinate data.
#' Only needed for non-spatial inputs.
#' @param y_coords String; name of column holding y-coordinate data. Only needed for non-spatial
#' inputs.
#' @param covariates RasterStack or RasterBrick; environmental covariate data.
#' @return A data frame with two attributes: soil class code and corresponding raster cell index
#' number.
#' @note The output of this function is an optional input for \code{\link{dsmartr_iterate}}. If
#' \code{known_points} is supplied as a data frame, the x and y coordinates must match the covariate
#' crs.
#' @examples \dontrun{
#' load('heronvale_known_sites')
#' load('heronvale_soilmap')
#' load('heronvale_covariates')
#'
#' # data frame input:
#' dsmartr_prep_points(src_map = heronvale_soilmap, known_points = heronvale_known_sites,
#' soil_id = 'CLASS', x_coords = 'x', y_coords = 'y', covariates = heronvale_covariates) }
#' @importFrom dplyr mutate_if
#' @importFrom raster cellFromXY crs extract
#' @importFrom sp spTransform
#' @importFrom sf st_as_sf st_set_crs st_transform
#' @export
dsmartr_prep_points <- function(known_points = NULL, soil_id  = NULL,
                                x_coords     = NULL, y_coords = NULL,
                                covariates   = NULL) {

  if (!dir.exists('inputs')) {
    dir.create('inputs', showWarnings = F)
  }
  in_dir <- file.path(getwd(), 'inputs')

  known_points <- mutate_if(known_points, is.factor, as.character)
  kpp <- if(is(known_points, 'sf')) {
    known_points <- st_transform(known_points, crs = crs(covariates, asText = TRUE))
    cellnos <- cellFromXY(covariates, xy = as(known_points, 'Spatial'))
    known_points$CELL <- cellnos
    known_points  <- st_set_geometry(known_points, NULL)
    names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
    known_points <- known_points[ , c('CLASS', 'CELL')]

    } else if (is(known_points, 'Spatial')) {
       known_points <- spTransform(known_points, CRSobj = crs(covariates))
       cellnos <- cellFromXY(covariates, xy = known_points)
       known_points <- known_points@data
       names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
       known_points$CELL  <- cellnos
       known_points <- known_points[ , c('CLASS', 'CELL')]

       } else if (is(known_points, 'data.frame')) {
         known_points <- st_as_sf(known_points, coords = c(x_coords, y_coords))
         known_points <- st_set_crs(known_points, crs = crs(covariates, asText = TRUE))
         cellnos   <- cellFromXY(covariates, xy = as(known_points, 'Spatial'))
         names(known_points)[names(known_points) == soil_id]  <- 'CLASS'
         known_points$CELL  <- cellnos
         known_points <- known_points[ , c('CLASS', 'CELL')]

       } else { stop('Please supply a valid input file') }

  saveRDS(kpp, file.path(in_dir, 'prepped_points.rds'))
  return(kpp)

  }
