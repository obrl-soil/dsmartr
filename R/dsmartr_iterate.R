#' Get all soil classes
#'
#' Returns the set of unique soil classes that occur in an input soil map, and, optionally,
#' a set of known points.
#' @keywords internal
#' @param soilmap Data Frame; Returned by \code{\link{prep_polygons}}.
#' @param soilpoints Data Frame; Returned by \code{\link{prep_points}}; optional.
#' @param col_stub String; stub of attribute column names holding soil class data
#' @return A factor holding the set of soil classes occurring in both datasets, sorted
#'   alphabetically. Internal to \code{\link{iterate}}.
#' @note It is up to the user to make sure the classification schema in both datasets matches/is
#'   compatible.
#' @examples \dontrun{
#' # run prep_polygons() and prep_points() examples, then:
#' class_levels <- get_classes(soil_map = pr_ap, soil_points = pr_points, col_stub = 'CLASS')}
#' @importFrom stats na.omit
get_classes <- function(soil_map = NULL, soil_points = NULL, col_stub = NULL) {

  map_levels <- unique(n_things(soil_map, col_stub))

  out_levels <- if(!is.null(soil_points)) {
    # sometimes known points have soil classes that have not been mapped
    point_levels <- unique(n_things(soil_points, col_stub))
    all_levels   <- base::union(map_levels, point_levels)
    as.factor(sort(all_levels))
  } else {
    as.factor(sort(map_levels))
  }
  out_levels
}

#' Sample a polygon for [iterate()]
#'
#' Randomly selects n cells for sampling and assigns them a soil class based on the overlying
#' map polygon's components
#' @keywords internal
#' @param pd A single-row sf object containing soil attributes and geometry. Usually output of
#'   running split(x, 1:nrow(x)) on a polygon sf dataframe.
#' @param cs String; stub of attribute column names holding soil class data
#' @param ps String; stub of attribute column names folding soil percentage data
#' @param nscol String; name of attribute holding number of samples needed
#' @param cellcol String; name of attriute holding list of raster cells to sample from
#' @param t_factor Integer; dirichlet dampener, supplied by parent function
#' @return A data frame containing two columns: weighted random allocation of soil classes,
#'   and cell numbers
#' @examples \dontrun{
#' # run prep_polygons() and prep_points() examples, then:
#' sample_points <- iter_sample_poly(pd = pr_ap[1, ], cs = 'CLASS', ps = 'PERC',
#' nscol = 'n_samples', cellcol = 'intersecting_cells', t_factor = t_factor)}
#' @importFrom gtools rdirichlet
#' @importFrom stats na.omit
iter_sample_poly <- function(pd = NULL, cs = NULL, ps = NULL,
                             nscol = NULL, cellcol = NULL, t_factor = NULL) {
  poly_nsample  <- unlist(pd[ , nscol], use.names = FALSE)
  poly_cells    <- unlist(pd[ , cellcol], use.names = FALSE)

  poly_cellsamp <- if(length(poly_cells) <= poly_nsample) {
    poly_cells
    } else {
      sample(poly_cells, size = poly_nsample, replace = FALSE)
    }

  poly_percs    <- as.numeric(n_things(pd, ps))
  poly_dirprops <- as.vector(rdirichlet(1, as.numeric(poly_percs) * t_factor))
  poly_classes  <- as.character(n_things(pd, cs))

  poly_alloc    <- mapply(function(class, dpn) {
      rep(class, times = dpn)
    },
    class = poly_classes,
    dpn   = ceiling(poly_dirprops * length(poly_cellsamp))
    )

  poly_alloc <- unlist(poly_alloc, use.names = FALSE)
  # shuffle randomly and make sure n is what it should be:
  # (sometimes you get n + 1 above)
  poly_alloc <- sample(poly_alloc, size = length(poly_cellsamp), replace = FALSE)
  # tried using a matrix here - speed boost insignificant and code harder to read
  poly_spoints  <- data.frame('CLASS' = poly_alloc,
                              'CELL'  = poly_cellsamp,
                              stringsAsFactors = FALSE)
}

#' Generate disaggregated soil maps
#'
#' Disaggregates an input soil map a given number of times.
#' @param prepped_map Data Frame; Returned by \code{\link{prep_polygons}}.
#' @param covariates RasterStack or RasterBrick; environmental covariate datasets.
#' @param prepped_points Data Frame; Returned by \code{\link{prep_points}}; optional.
#' @param id_field String; name of unique polygon identifier field.
#' @param n_iterations Integer; desired number of model runs.
#' @param t_factor Integer; dirichlet distribution modifier. Where map units have >1 soil class
#'   component, larger values of \code{t_factor} dampen variation away from the input proportions.
#' @param c5_ctrl List; output of \code{\link[C50]{C5.0Control}}; optional.
#' @param cpus Integer; number of processors to use in parallel.
#' @param write_files String; choose 'native' to write outputs in GeoTIFF, GPKG, or csv as
#' appropriate, or 'rds' for outputs in RDS format. Note that map rasters are always
#' written as GeoTIFF.
#' @param write_samples Logical; whether to retain the covariate samples taken during each
#'   iteration. If \code{TRUE}, each set of samples is written to disk as GPKG or RDS.
#' @param resume_from Integer; If this function is interrupted, it can be resumed from a specified
#'   iteration number (prevents overwriting existing model iterations in output directory).
#' @note This function writes a large number of files to disk.
#' @return RasterStack; \code{n_iterations} layers. As a side effect, for each model run, outputs
#'   written to disk include the predicted soils map, the C5.0 model summary, and optionally, the
#'   covariate sample points in spatial data format.
#' @examples \dontrun{
#' # Polygons only:
#' # run prep_polygons() example code, then
#' iteration_maps <- iterate(prepped_map = pr_ap, covariates = heronvale_covariates,
#'  id_field = 'POLY_NO', n_iterations = 20,
#'  cpus = max(1, (parallel::detectCores() - 1)), write_files = 'native', write_samples = TRUE)
#'
#'  # Polygons, points and a C50 model tweak, no samples, rds output:
#'  win_on <- C50::C5.0Control(winnow = TRUE)
#'  iteration_maps <- iterate(prepped_map = prepped_ap, covariates = heronvale_covariates,
#'  id_field = 'POLY_NO', prepped_points = pr_pts, n_iterations = 20, c5_ctrl = win_on,
#'  cpus = max(1, (parallel::detectCores() - 1)), write_files = 'rds')
#'
#' ## Oh no, there was a blackout halfway through my model run:
#' iteration_maps_2 <- iterate(prepped_map = pr_ap, covariates = heronvale_covariates,
#' id_field = 'POLY_NO', n_iterations = 6, cpus = max(1, (parallel::detectCores() - 1)),
#' write_files = 'all', write_samples = TRUE, resume_from = 14)
#' }
#' @importFrom C50 C5.0
#' @importFrom dplyr distinct filter
#' @importFrom gtools rdirichlet
#' @importFrom purrr map
#' @importFrom raster beginCluster clusterR endCluster extract inMemory readAll writeRaster xyFromCell
#' @importFrom sf st_as_sf st_point st_sfc st_write
#' @importFrom stats complete.cases na.omit
#' @importFrom tidyr gather
#' @importFrom utils capture.output setTxtProgressBar txtProgressBar write.table
#' @export
iterate <- function(prepped_map    = NULL,
                    covariates     = NULL,
                    prepped_points = NULL,
                    id_field       = NULL,
                    n_iterations   = NULL,
                    t_factor       = 10000,
                    c5_ctrl        = NULL,
                    cpus           = 1,
                    write_files    = c('native', 'rds'),
                    write_samples  = FALSE,
                    resume_from    = NULL) {

  # output directories
  dir.create(file.path(getwd(), 'iterations', 'maps'),   showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(getwd(), 'iterations', 'models'), showWarnings = FALSE, recursive = TRUE)
  strr   <- file.path(getwd(), 'iterations', 'maps')
  strm   <- file.path(getwd(), 'iterations', 'models')

  class_levels <- get_classes(soil_map    = prepped_map,
                              soil_points = prepped_points,
                              col_stub   = 'CLASS')

  message(paste0(Sys.time(), ': dsmartr iteration in progress...'))
  pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

  ## realisations
  start <- if (is.null(resume_from)) { 1 } else { resume_from }

  lapply(seq.int(start, n_iterations), function(j) {

    beginCluster(cpus)

    src_split <- split(prepped_map, 1:nrow(prepped_map))
    sample_points <- map(.x = src_split, .f = function(z) {
      iter_sample_poly(pd = z, cs = 'CLASS', ps = 'PERC', nscol = 'n_samples',
                       cellcol = 'intersecting_cells', t_factor = t_factor)})
    all_samplepoints <- do.call('rbind', sample_points)

    # sample covariates
    cov_sample <- raster::extract(covariates, y = all_samplepoints$CELL)  ## BOTTLENECK
    all_samplepoints <- cbind(all_samplepoints, cov_sample)

    # optionally, add known locations
    all_samplepoints <- if(!is.null(prepped_points)) {

      known_sample   <- raster::extract(covariates, y = prepped_points$CELL)
      prepped_points <- cbind(prepped_points, known_sample)

      # make sure known cells haven't been randomly sampled in this iteration
      # if so, drop the randomly generated row in favour of the known class
      asp <- dplyr::filter(all_samplepoints, !(CELL %in% c(prepped_points$CELL)))
      asp <- rbind(asp, prepped_points)
      asp
    } else {
      all_samplepoints
    }

    # force all outputs to be on the same scale eg. output map value 1 always equals factor level 1
    all_samplepoints$CLASS         <- as.factor(all_samplepoints$CLASS)
    levels(all_samplepoints$CLASS) <- class_levels

    # set up model input
    model_input <- all_samplepoints[complete.cases(all_samplepoints), ]

    # generate decision tree
    res <- if(is.null(c5_ctrl)) {
      C5.0(x = model_input[, !(names(model_input) %in% c('CLASS', 'CELL'))],
           y = model_input$CLASS)
    } else {
      C5.0(x = model_input[, !(names(model_input) %in% c('CLASS', 'CELL'))],
           y = model_input$CLASS,
           control = c5_ctrl)
    }

    # make prediction map
    iter_map <- clusterR(covariates, raster::predict, args = list(res), datatype = 'INT2S')

    # factorise prediction map
    lookup <- data.frame("ID"    = as.integer(as.factor(res$levels)),
                         "CLASS" = as.factor(res$levels))
    levels(iter_map) <- lookup

    # write probability map from this realisation to GeoTIFF
    # (lookup values are in accompanying aux file)
    writeRaster(iter_map,
                filename  = file.path(strr, paste0('map_', j, '.tif')),
                format    = "GTiff",
                datatype  = "INT2S",
                NAflag    = -9999,
                overwrite = TRUE)

    # may as well write lookup to csv as aux.xml files aren't read by QGIS :(
    if (j == start) {
      # NB Excel is the default csv viewer for many people but it doesn't like csv files
      # where column 1 is called 'ID'
      names(lookup) <- c('VALUE', 'CLASS')
      write.table(lookup, file.path(strr, 'class_lookup.csv'),
                  col.names = TRUE,  row.names = FALSE,
                  quote     = FALSE, sep       = ",")
    }

    # write model from this realisation to file (summary only in txt)
    if (write_files == 'rds') {
      saveRDS(res, file.path(strm, paste0('C5_model_', j, '.rds')))
    } else {
      out <- capture.output(summary(res))
      f2  <- file.path(strm, paste0('C5_model_', j, '.txt'))
      cat(out, file = f2, sep = "\n", append = TRUE)
      }

    if (write_samples == TRUE) {

      if (!dir.exists(file.path(getwd(), 'iterations', 'samples'))) {
        dir.create(file.path(getwd(), 'iterations', 'samples'),
                   showWarnings = FALSE, recursive = TRUE)
      }
      strd <- file.path(getwd(), 'iterations', 'samples')

      # spatialise (sf point object)
      model_input$geometry <- st_sfc(lapply(model_input$CELL, function(x) {
        st_point(as.vector(xyFromCell(covariates, x)))
      }))
      model_input <- st_as_sf(model_input,
                              crs = covariates@crs@projargs, agr = 'constant')

      if (write_files == 'rds') {
        saveRDS(model_input, file.path(strd, paste0('samples_',  j, '.rds')))
      } else {
        suppressWarnings(st_write(obj          = model_input,
                                  dsn          = file.path(strd, paste0('samples_', j, '.gpkg')),
                                  driver       = 'GPKG',
                                  delete_dsn   = TRUE,
                                  delete_layer = TRUE,
                                  quiet        = TRUE)
                         )
      }
    }

    endCluster()
    setTxtProgressBar(pb, j)
  } )

  # return maps to for dsmartr_collate
  map_list  <- list.files(strr, pattern = '\\.tif$', full.names = TRUE)
  iteration_maps <- raster::stack(map_list)
  names(iteration_maps) <- paste0('iter_', 1:length(map_list))


  close(pb)
  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'iterations')))

  iteration_maps

}
