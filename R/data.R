#' Demonstration covariates for dsmartr
#'
#' A dataset containing a variety of environmental parameters that may correlate with mapped soil
#' concepts.
#' @format A RasterBrick with dimensions 299, 285, 85215, 13  (nrow, ncol, ncell, nlayers) and
#' resolution of 26.6 m, projected in Albers Equal Area (EPSG:3577). The dataset is clipped to
#' the extent of the \href{https://data.qld.gov.au/dataset/queensland-map-sheet-key-maps-series}{Heronvale 1:10000 topographic key map sheet}, located on the coast of Northern Queensland
#' between Proserpine and Bowen. Included surfaces are:
#' \describe{
#'   \item{DEM_S}{GeoScience Australia's
#'   \href{http://www.ga.gov.au/metadata-gateway/metadata/record/72759/}{1" SRTM DEM-S } product}
#'   \item{GM_INT_pattern}{'Intensity' output of GRASS add-on module
#'   \href{https://grass.osgeo.org/grass72/manuals/addons/r.geomorphon.html}{r.geomorphon}, calculated
#'   from DEM_S in GRASS 7.2 with a search distance of 11 cells, a skip distance of 3 cells, and a
#'   flatness threshold of 0.3 degrees. This roughly corresponds to the 'landscape pattern' concept
#'   outlined in \href{http://www.publish.csiro.au/book/5230}{NCST, 2009}. 'Intensity' is
#'   the mean difference in elevation between a pixel and the eight neighbourhood pixels used to
#'   define its geomorphon.}
#'   \item{GM_MCF_pattern}{'Most common forms' (MCF) output of r.geomorphon, settings as above. MCF
#'   is a 10-category set of simple landscape shapes - flat, summit, ridge, shoulder, spur, slope,
#'    hollow, footslope, valley, depression.}
#'    \item{GM_RNG_pattern}{'Range' output of r.geomorphon, settings as above. 'Range' is the
#'    difference between the max and min elevations of the cells used to define the geomorphon.}
#'    \item{NETD_MRVBF}{Multi-Resolution Valley Bottom Flatness Index derived from 1" DEM-S,
#'    published by the CSIRO and available from \url{http://doi.org/10.4225/08/5701C885AB4FE}}
#'    \item{NETD_Prescott}{Prescott Index derived from 1" DEM-S, published by the CSIRO and
#'    available from \url{http://doi.org/10.4225/08/53EB2D0EAE377}}
#'    \item{NETD_SLOPE_DG}{Slope derived from 1" DEM-S, published by the CSIRO and available from
#'    \url{http://doi.org/10.4225/08/5689DA774564A}}
#'    \item{PCRE_DBVG1M}{\href{https://data.qld.gov.au/dataset/pre-clearing-broad-vegetation-groups-of-queensland-version-3-0}{Pre-clearing
#'    dominant broad vegetation groups of Queensland, v3.0}. This dataset was
#'    rasterised from a polygon source.}
#'    \item{RM_F_K}{\href{http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/rr2/National_Coverages/http/radmap_v3_2015_filtered_pctk/catalog.xml}{Radiometric
#'    Map of Australia v3 - Filtered percent Potassium}}
#'    \item{RM_F_T}{\href{http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/rr2/National_Coverages/http/radmap_v3_2015_filtered_ppmth/catalog.xml}{Radiometric
#'    Map of Australia v3 - Filtered ppm Thorium}}
#'    \item{RM_F_U}{\href{http://dap.nci.org.au/thredds/remoteCatalogService?catalog=http://dapds00.nci.org.au/thredds/catalog/rr2/National_Coverages/http/radmap_v3_2015_filtered_ppmu/catalog.xml}{Radiometric
#'    Map of Australia v3 - Filtered ppm Uranium}}
#'    \item{SAGA_MSPOS}{\href{http://www.saga-gis.org/saga_tool_doc/2.3.0/ta_morphometry_14.html}{SAGA
#'    Mid-slope position}, calculated from DEM-S in SAGA 2.3.}
#'    \item{SAGA_WI}{\href{http://www.saga-gis.org/saga_tool_doc/2.3.0/ta_hydrology_15.html}{SAGA
#'    Wetness Index}, calculated from DEM-H in SAGA 2.3.}}
#' All covariates were processed in GRASS 7 to align with the DEM-S, reprojected to EPSG:3577 with
#' appropriate resampling (bilinear for continuous data, nearest-neighbour for categorical).
'heronvale_covariates'

#' Points where soil class is known
#'
#' A dataset containing the locations of described soil profiles within the Heronvale 1:10000
#' topographic key map sheet. Sites are from a soil and land resource survey conducted in the early
#' 2000s. Full profile descriptions can be viewed in the
#' \href{https://qldglobe.information.qld.gov.au/}{Queensland Globe}.
#' @format An tbl/sf_POINT object in EPSG:3577 with 36 rows and 6 columns:
#' \describe{
#'  \item{PROJECT_CODE}{Three-letter code identifying the soil survey}
#'  \item{SITE_ID}{Within-project unique site identifier}
#'  \item{LONGITUDE}{Longitude coordinate in EPSG:4283 (GDA94)}
#'  \item{LATITUDE}{Latitude coordinate in EPSG:4283 (GDA94)}
#'  \item{CLASS}{Soil class identified at this site}
#'  \item{geom}{list-column with class 'sfc_GEOMETRY', 'sfc'}}
#'  Site locations are accurate to within +/-50m.
#' @source \itemize{
#' \item{\url{https://publications.qld.gov.au/dataset/soils-whitsunday-coast-wcs}}
#' \item{\url{https://www.qld.gov.au/environment/land/soil/soil-data/survey-types/}}}
'heronvale_known_sites'

#' Broadscale soil map of the Heronvale area
#'
#' A dataset containing an extract of 1:100,000 scale soils mapping from the Whitsunday coast area.
#' The complete dataset can be viewed in the \href{https://qldglobe.information.qld.gov.au/}{Queensland
#'  Globe.} Data was extracted from the Queensland Soil and Land Information (SALI) database in 2017
#'  and formatted for dsmartr in R 3.3.
#'  @format An tbl/sf dataframe of mixed POLYGON/MULTIPOLYGON types in EPSG:3577, with 105 rows and
#'  10 variables:
#'  \describe{
#'  \item{PROJECT_CODE}{Three-letter code identifying the soil survey}
#'  \item{POLY_NO}{Within-project Unique Map Area identifier}
#'  \item{POLY_UID}{Concatenation of PROJECT_CODE and POLY_NO}
#'  \item{CLASS_1}{dominant soil class within the map unit}
#'  \item{CLASS_2}{subdominant soil class within the map unit}
#'  \item{CLASS_3}{subdominant soil class within the map unit}
#'  \item{PERC_1}{Estimated percentage of the map unit occupied by CLASS_1}
#'  \item{PERC_2}{Estimated percentage of the map unit occupied by CLASS_2}
#'  \item{PERC_3}{Estimated percentage of the map unit occupied by CLASS_3}
#'  \item{geom}{list-column with class 'sfc_GEOMETRY', 'sfc'}}
#' @source \itemize{
#' \item{\url{https://publications.qld.gov.au/dataset/soils-whitsunday-coast-wcs}}
#' \item{\url{https://www.qld.gov.au/environment/land/soil/soil-data/survey-types/}}}
'heronvale_soilmap'

#' Decoded soil class names
#'
#' A table of soil class labels found in 'heronvale_soilmap' and 'heronvale_known_sites', along with
#' their full names. Explanations of the concept associated with each soil name are available in the
#'  report associated with this mapping project.
#' @format A tibble with 38 rows and 2 columns:
#' \describe{
#' \item{CLASS}{Two-letter soil class label}
#' \item{NAME}{Soil class name}}
#' @source \itemize{
#' \item{\url{https://publications.qld.gov.au/dataset/soils-whitsunday-coast-wcs}}
#' \item{\url{https://www.qld.gov.au/environment/land/soil/soil-data/survey-types/}}}
'heronvale_soilnames'
