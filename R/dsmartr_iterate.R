#' dsmartr model iterator
#'
#' Disaggregates an input soil map a given number of times.
#' @param prepped_map Data Frame; Returned by \code{\link{dsmartr_prep_polygons}}.
#' @param covariates RasterStack or RasterBrick; environmental covariate datasets.
#' @param prepped_points Data Frame; Returned by \code{\link{dsmartr_prep_points}}; optional.
#' @param id_field String; name of unique polygon identifier field.
#' @param n_iterations Integer; desired number of model runs.
#' @param t_factor Integer; dirichlet distribution modifier. Where map units have >1 soil class
#' component, larger values of \code{t_factor} dampen variation away from the input proportions.
#' @param c5_ctrl List; output of \code{\link[C50]{C5.0Control}}; optional.
#' @param cpus Integer; number of processors to use in parallel.
#' @param write_files String; choose 'all' to write outputs in both rds and 'normal' files (GeoTIFF,
#'  ESRI shapefile, csv, txt as appropriate).
#' @param write_samples Logical; whether to retain the covariate samples taken during each
#' iteration. If \code{TRUE}, each set of samples is written to disk as rds and/or ESRI Shapefile.
#' @param resume_from Integer; If this function is interrupted, it can be resumed from a specified
#' iteration number (prevents overwriting existing model iterations in output directory).
#' @note This function writes a large number of files to disk.
#' @return RasterStack; \code{n_iterations} layers. As a side effect, for each model run, outputs
#' written to disk include the predicted soils map, the C5 model summary, and optionally, the
#' covariate sample points in spatial data format.
#' @examples \dontrun{
#' ## Polygons only:
#' # run dsmartr_prep_polygons() and then
#' dsmartr_iterate(prepped_map = prepped_polygons, covariates = heronvale_covariates,
#'  id_field = 'POLY_NO', n_iterations = 20,
#'  cpus = max(1, (parallel::detectCores() - 1)), write_files = 'all', write_samples = TRUE)
#'
#' ## Oh no, there was a blackout halfway through my model run:
#' dsmartr_iterate(prepped_map = prepped_polygons, covariates = heronvale_covariates,
#' id_field = 'POLY_NO', n_iterations = 6, cpus = max(1, (parallel::detectCores() - 1)),
#' write_files = 'all', write_samples = TRUE, resume_from = 14)}
#' @importFrom tidyr gather
#' @importFrom dplyr distinct filter
#' @importFrom purrr map
#' @importFrom gtools rdirichlet
#' @importFrom raster extract
#' @importFrom C50 C5.0
#' @export
dsmartr_iterate <- function(prepped_map    = NULL,
                            covariates     = NULL,
                            prepped_points = NULL,
                            id_field       = NULL,
                            n_iterations   = NULL,
                            t_factor       = 10000,
                            c5_ctrl        = NULL,
                            cpus           = 1,
                            write_files    = c('rds_only', 'all'),
                            write_samples  = FALSE,
                            resume_from    = NULL) {

  # output directories
  dir.create('iterations',        showWarnings = F)
  dir.create('iterations/maps',   showWarnings = F)
  dir.create('iterations/models', showWarnings = F)
  strr   <- file.path(getwd(), 'iterations', 'maps')
  strm   <- file.path(getwd(), 'iterations', 'models')

  # set consistent factoring across model runs
  prepped_map <- mutate_if(prepped_map, is.factor, as.character)
  if(!is.null(prepped_points)) {
    # sometimes known points have soil classes that have not been mapped
    class_levels <- prepped_map[, c(grep('CLASS_', names(prepped_map)))]
    class_levels <- gather(data = class_levels, na.rm = TRUE)
    class_levels <- na.omit(distinct(class_levels, value))
    class_levels <- unlist(class_levels, use.names = FALSE)
    point_levels <- as.character(unique(unlist(prepped_points$CLASS)))
    class_levels <- as.factor(sort(union(class_levels, point_levels)))
  } else {
    class_levels <- prepped_map[, c(grep('CLASS_', names(prepped_map)))]
    class_levels <- gather(data = class_levels, na.rm = TRUE)
    class_levels <- na.omit(distinct(class_levels, value))
    class_levels <- unlist(class_levels, use.names = FALSE)
    class_levels <- as.factor(sort(class_levels))
  }

  message(paste0(Sys.time(), ': dsmartr iteration in progress...'))
  pb <- txtProgressBar(min = 0, max = n_iterations, style = 3)

  ## realisations
  start <- if (is.null(resume_from)) { 1 } else { resume_from }

  lapply(seq.int(start, n_iterations), function(j) {

    beginCluster(cpus)

    src_split <- split(prepped_map, 1:nrow(prepped_map))
    sample_points <- map(.x = src_split, .f = function(z) {
        poly_cells    <- unlist(z[ , 'intersecting_cells'], use.names = FALSE)
        poly_nsample  <- unlist(z[ , 'n_samples'], use.names = FALSE)
        poly_cellsamp <- if (length(poly_cells) <= poly_nsample) {
          poly_cells
        } else {
          sample(poly_cells, size = poly_nsample, replace = FALSE)
        }
        poly_percs    <- as.vector(na.omit(unlist(z[, c(grep('PERC_', names(z)))])))
        poly_dirprops <- as.vector(rdirichlet(1, as.numeric(poly_percs) * t_factor))
        poly_classes  <- as.vector(na.omit(unlist(z[, c(grep('CLASS_', names(z)))])))
        poly_alloc    <- mapply(function(class, dpn) {
          rep(class, times = dpn)
        },
        class = poly_classes,
        dpn   = ceiling(poly_dirprops * length(poly_cellsamp))
        )
        poly_alloc <- unlist(poly_alloc)
        # shuffle randomly and make sure n is what it should be:
        #(sometimes you get n + 1 above)
        poly_alloc <- sample(poly_alloc, size = length(poly_cellsamp), replace = FALSE)
        # tried using a matrix here - speed boost insignificant and code harder to read
        poly_spoints  <- data.frame('CLASS' = poly_alloc,
                                    'CELL'  = poly_cellsamp)
      })

    # get all the sampling data for all the polygons for this iteration into one object
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
    asp <- filter(all_samplepoints, !(CELL %in% c(prepped_points$CELL)))
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

    # factorise iter_map
    lookup <- data.frame("ID"    = as.integer(as.factor(res$levels)),
                         "CLASS" = as.factor(res$levels))
    levels(iter_map) <- lookup

    if (write_files == 'rds_only') {

      saveRDS(res, file.path(strm, paste0('C5_model_', j, '.rds')))
      # stop file reference pointing at tempdir gri/grd as they don't persist
      mout <- if(inMemory(iter_map) == FALSE) { readAll(iter_map) } else { iter_map }
      mout@file@name <- file.path(strr, paste0('map_', j, '.rds'))
      suppressWarnings(saveRDS(mout, file.path(strr, paste0('map_', j, '.rds'))))

    } else {

      saveRDS(res, file.path(strm, paste0('C5_model_', j, '.rds')))
      mout <- if(inMemory(iter_map) == FALSE) { readAll(iter_map) } else { iter_map }
      mout@file@name <- file.path(strr, paste0('map_', j, '.rds'))
      suppressWarnings(saveRDS(mout, file.path(strr, paste0('map_', j, '.rds'))))

      # write decision tree from this realisation to text
      out <- capture.output(summary(res))
      f2  <- file.path(strm, paste0('C5_model_', j, '.txt'))
      cat(out, file = f2, sep = "\n", append = TRUE)

      # write probability map from this realisation to GeoTIFF (lookup values are embedded)
      nme <- file.path(strr, paste0('map_', j, '.tif'))
      writeRaster(iter_map, filename = nme, format = "GTiff", overwrite = T, datatype = "INT2S",
                  NAflag = -9999)

      # may as well write lookup to csv as aux.xml files aren't read by QGIS :(
      if (j == 1) {
        # NB Excel is the default csv viewer for many people but it doesn't like csv files
        # where column 1 is called 'ID'
        names(lookup) <- c('VALUE', 'CLASS')
        write.table(lookup, file.path(strr, 'class_lookup.csv'),
                    col.names = TRUE,  row.names = FALSE,
                    quote     = FALSE, sep       = ",")
      }
    }

    if (write_samples == TRUE) {
      dir.create('iterations/samples', showWarnings = F)
      strd <- file.path(getwd(), 'iterations', 'samples')
      # spatialise (sf point object)
      allsamp_sf <- all_samplepoints
      allsamp_sf$geometry <- st_sfc(lapply(allsamp_sf$CELL, function(x) {
        st_point(as.vector(xyFromCell(covariates, x)))
        }))
      allsamp_sf <- st_as_sf(allsamp_sf, crs = covariates@crs@projargs, agr = 'contstant')

      if (write_files == 'rds_only') {
        saveRDS(allsamp_sf, file.path(strd, paste0('samples_',  j, '.rds')))

      } else {

        saveRDS(allsamp_sf, file.path(strd, paste0('samples_',  j, '.rds')))

        ### NB Don't change this to GPKG format yet - disk write time is crazy slow

        # make a lookup table for all_samplepoints covariate column names, because they're about
        # to get severely abbreviated for writing to shp
        cov_LUT_nm   <- file.path(strd, 'covariate_LUT.csv')
        cov_names    <- names(all_samplepoints[5:ncol(all_samplepoints)])
        cov_shpnames <- paste0('COV_', 1:length(cov_names))
        if (!file.exists(cov_LUT_nm)) {
          cov_LUT <- data.frame("COV_NAMES" = cov_names, "SHPCOL_NAMES" = cov_shpnames)
          write.table(cov_LUT, file = cov_LUT_nm, sep = ', ',
                      quote = FALSE, col.names = TRUE, row.names = FALSE) }

        names(allsamp_sf)[5:(ncol(allsamp_sf) - 1)] <- cov_shpnames
        sf_name    <- paste0('samples_', j, '.shp')
        # see https://trac.osgeo.org/gdal/ticket/6803,
        # and https://github.com/edzer/sfr/issues/306:
        suppressWarnings(write_sf(allsamp_sf,
                 file.path(strd, sf_name),
                 driver = 'ESRI Shapefile',
                 delete_layer = TRUE))
      }
    }

    endCluster()
    setTxtProgressBar(pb, j)
  } )

  # send maps to global env for dsmartr_collate
  map_rds  <- list.files(strr, pattern = '\\.rds$', full.names = TRUE)
  map_list <- lapply(map_rds, function(x) readRDS(x))
  names(map_list) <- paste0('iter_', 1:length(map_list)) # still a list
  assign('iteration_maps', raster::stack(map_list), envir = .GlobalEnv)

  close(pb)
  message(paste0(Sys.time(), ': ...complete. dsmartr outputs can be located at ',
                 file.path(getwd(), 'iterations')))

}
