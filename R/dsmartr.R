#' dsmartr: An R implementation of the DSMART algorithm
#'
#' dsmartr is an approach to the DSMART model (Odgers, 2014), which attempts increase the
#' level of detail in legacy soil maps. Relationships between these maps and various environmental
#' datasets are used to build a classification tree model. The model is used to map the map
#' components as a fine-scaled raster. This process is repeated multiple times, with variation introduced around sample location and
#' soil class proportions to vary the outputs somewhat. The model realisations are then processed
#' to give a 'most likely' soils map. Refer to the package vignette for more detail.
#'
#' @references Odgers, N.P., Sun, W., McBratney, A.B., Minasny, B., Clifford, D., (2014)
#' \href{http://dx.doi.org/10.1016/j.geoderma.2013.09.024}{Disaggregating and harmonising soil map
#' units through resampled classification trees}. Geoderma, 214-215: 91-100.
#'
"_PACKAGE"
