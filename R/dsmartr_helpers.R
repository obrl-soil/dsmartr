#' helper function
#'
#' Returns cell values from selected columns in an attribute table as an atomic vector,
#' with NA values excluded.
#' @keywords internal
#' @param input An object of class sfc_POLYGON/MULTIPOLYGON, see
#' \code{\link{dsmartr_prep_polygons}}.
#' @param selector String; a column name or name stub.
#' @return An atomic vector of non-NA values in the selected columns, by order of occurrence.
#' @note This is a helper function and not widely applicable. It is expected to be
#'   used by-row on wide-formatted sf attribute data.
#' @examples \dontrun{
#' load('heronvale_soilmap')
#' percs_1 <- n_things(heronvale_soilmap[1, ], 'PERC')}
#' @importFrom stats na.omit
n_things <- function(input = NULL, selector = NULL) {
  input <- as.data.frame(input, stringsAsFactors = FALSE)
  if(length(grep(selector, names(input))) == 0) {
    stop("Error: Selector not found within input's column names")
  } else {
  output <- as.vector(na.omit(unlist(input[, c(grep(selector, names(input)))])))
  if(length(output) == 0) { NA } else { output }}
}


#' Check dsmartr polygon inputs
#'
#' Checks dsmartr input polygons for potential data quality issues.
#' @param src_map An sfc_POLYGON/MULTIPOLYGON or SpatialPolygonsDataFrame object.
#' @param id_field String; name of unique identifier field.
#' @param cs String; stub of attribute column names holding soil class data
#' @param ps String; stub of attribute column names folding soil percentage data
#' @details This function highlights some issues that may cause dsmartr to either fail or behave
#'   unpredictably. Usually these are side-effects of prepration work carried out in other software.
#'   Sometimes they are the result of errors in the source data. Note that polygon geometry is not
#'   explicitly checked, only attributes. At present, four problems are checked - whether data is
#'   missing (e.g. a polygon class without a matching percentage), whether a percentage is 0%,
#'   whether total percentage (by polygon) is not 100%, and whether any IDs are duplicated.
#' @return An sf object with appended logical columns indicating presence of possible data faults.
#' @examples \dontrun{
#' load("heronvale_soilmap")
#' checked_map <- dsmartr_check_polygons(src_map = heronvale_soilmap,
#' id_field = 'POLY_NO')}
#' @export
dsmartr_check_polygons <- function(src_map  = NULL, id_field = NULL,
                                   cs = NULL, ps = NULL) {

  # coerce to sf
  src_map <- if(is(src_map, 'sf') == FALSE) {
    st_as_sf(src_map) # numpty-proof this
  } else {
    src_map
  }

  ### highlight problems in new columns for easy ID
  src_split <- split(src_map, 1:nrow(src_map))

  missing_data <- purrr::map_lgl(src_split, function(md) {
               n_classes <- length(n_things(md, cs))
               n_percs   <- length(n_things(md, ps))
               ifelse(n_classes != n_percs, TRUE, FALSE)
               })

  zero_percs <- purrr::map_lgl(src_split, function(zp) {
               poly_percs <- n_things(zp, ps)
               ifelse(any(poly_percs == 0), TRUE, FALSE)
             })

  problem_percs <- purrr::map_lgl(src_split, function(sps) {
               total <- sum(n_things(sps, ps))
               ifelse(total != 100, TRUE, FALSE)
               })

  src_map$missing_data  <- as.vector(missing_data)
  src_map$zero_percs    <- as.vector(zero_percs)
  src_map$problem_percs <- as.vector(problem_percs)
  src_map$duplicate_ids <- duplicated(src_map[, id_field])

  return(src_map)
}

#' generate masks from samples
#'
#' Generates a rasterlayer for each dsmart iteration highlighting areas where
#' one or more covariates has a value that is out of range of the sample used to
#' model soil distribution.
#' @param samples List; POINT_sfc objects created by \code{\link{dsmartr_iterate}} when option
#' write_samples = TRUE
#' @param covariates RasterStack or RasterBrick; environmental covariate datasets.
#' @param tolerance Integer; the number of out of range covariates allowed. Defaults to 0,
#' that is, if any covariate is out of range on a given pixel, that pixel will be masked.
#' @param cpus Integer; number of processors to use in parallel.
#' @details For each model iteration this function generates a raster that can be used to mask
#'  dsmartr predictions where covariate values are out of range of the sample of cells
#'  used in that iteration. The masks can also be combined to produce an overall mask, applicable to
#'  final products.
#' @return a list of rasters
#' @examples \dontrun{
#' # run dsmartr_iterate() with the example data, then:
#' sample_list <- lapply(as.list(list.files(file.path(getwd(), 'iterations', 'samples'),
#' pattern = '\\.shp', full.names = TRUE)), function(x) read_sf(x))
#' masks <- dsmartr_pred_masks(samples = sample_list,
#' covariates = heronvale_covariates, cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # less strict: up to 3 out of range covariates are allowed
#' masks <- dsmartr_pred_masks(samples = sample_list, covariates = heronvale_covariates,
#' tolerance = 3, cpus = max(1, (parallel::detectCores() - 1)))
#'
#' # a 'final product' mask, where outputs are masked only if the mode of that
#' cell is 'do not predict'.
#' all_masks <- stack(masks)
#' modal_mask <- calc(all_masks, function(cell) ifelse(modal(cell) == FALSE, NA ,1))
#' # the above could be applied to the most-likely soils map after running dsmartr_most_likely:
#' masked_m1 <- most_likely[[1]][[1]] + modal_mask
#' }
#' @importFrom raster calc writeRaster
#' @importFrom sf st_set_geometry
#' @importFrom purrr map2_lgl
#' @export
dsmartr_pred_masks <- function(samples = NULL, covariates = NULL, tolerance = 0L, cpus = 1) {
  if (!dir.exists(file.path(getwd(), 'iterations', 'masks'))) {
    dir.create(file.path(getwd(), 'iterations', 'masks'),
               showWarnings = FALSE, recursive = TRUE)
  }
  drmsk <- file.path(getwd(), 'iterations', 'masks')

  message(paste0(Sys.time(), ': dsmartr prediction mask creation in progress...'))
  pb <- txtProgressBar(min = 0, max = length(samples), style = 3)

  all_masks <- mapply(FUN = function(m, n) {
    sx <- st_set_geometry(m, NULL)
    rangex <- apply(sx[, c(3:ncol(sx))], MARGIN = 2,
                    FUN = function(x) range(x, na.rm = TRUE))
    beginCluster(n = cpus)
    pred_mask <- clusterR(x    = covariates,
                          fun  = calc,
                          args = list(function(cell) {
                            in_or_out <- map2_lgl(.x = cell, .y = seq_along(cell),
                                                  function(.x, .y) {
                                                    .x >= rangex[1, .y] & .x <= rangex[2, .y]
                                                    })
                            ifelse(sum(in_or_out == FALSE) > tolerance, NA, 0)
                            }),
                          filename = file.path(drmsk, paste0('pred_mask_', n, '.tif')),
                          NAflag = -9999,
                          datatype = 'INT2S',
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
