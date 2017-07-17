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
  as.vector(na.omit(unlist(as.data.frame(input)[, c(grep(selector, names(input)))])))
}

#' Check dsmartr polygon inputs
#'
#' Checks dsmartr input polygons for potential data quality issues.
#' @param src_map An sfc_POLYGON/MULTIPOLYGON or SpatialPolygonsDataFrame object.
#' @param id_field String; name of unique identifier field.
#' @details This function highlights some issues that may cause dsmartr to either fail or behave
#'   unpredictably. Usually these are side-effects of prepration work carried out in other software.
#'   Sometimes they are the result of errors in the source data. Note that polygon geometry is not
#'   explicitly checked, only attributes. At present, three problems are checked - whether data is
#'   missing (e.g. a polygon class without a matching percentage), whether a percentage is 0%, and
#'   whether total percentage (by polygon) is not 100%.
#' @return An sf object with appended logical columns indicating presence of possible data faults.
#' @examples \dontrun{
#' load("heronvale_soilmap")
#' checked_map <- dsmartr_check_polygons(src_map = heronvale_soilmap,
#' id_field = 'POLY_NO')}
#' @export
dsmartr_check_polygons <- function(src_map  = NULL, id_field = NULL) {

  # coerce to sf
  src_map <- if(is(src_map, 'sf') == FALSE) {
    st_as_sf(src_map) # numpty-proof this
  } else {
    src_map
  }

  ### highlight problems in new columns for easy ID
  src_split <- split(src_map, 1:nrow(src_map))

  missing_data <- purrr::map_lgl(src_split, function(md) {
               n_classes <- length(n_things(md, 'CLASS_'))
               n_percs   <- length(n_things(md, 'PERC_'))
               ifelse(n_classes != n_percs, TRUE, FALSE)
               })

  zero_percs <- purrr::map_lgl(src_split, function(zp) {
               poly_percs <- n_things(zp, 'PERC_')
               ifelse(any(poly_percs == 0), TRUE, FALSE)
             })

  problem_percs <- purrr::map_lgl(src_split, function(sps) {
               total <- sum(n_things(sps, 'PERC_'))
               ifelse(total != 100, TRUE, FALSE)
               })

  src_map$missing_data  <- as.vector(missing_data)
  src_map$zero_percs    <- as.vector(zero_percs)
  src_map$problem_percs <- as.vector(problem_percs)

  return(src_map)
}
