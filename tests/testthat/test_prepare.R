context('dsmartr_prepare')

data('heronvale_soilmap')
data('heronvale_covariates')

# test on sf object
test_that('prep_polygons with flat rate', {
          c(
            test_1 <- prep_polygons(src_map       = heronvale_soilmap[c(1L,2L,11L,75L, 83L), ],
                                    covariates    = heronvale_covariates,
                                    id_field      = 'POLY_UID',
                                    sample_method = 'flat',
                                    sample_rate     = 10L),
            expect_output(str(test_1), '5 obs'),
            expect_output(str(test_1), '13 variables'),
            expect_s3_class(test_1$area_sqkm, 'units'),
            expect_equal(test_1$n_soils, c(2,3,1,3,2)),
            expect_equal(test_1$n_samples, c(10, 10, 10, 10, 10)),
            expect_is(test_1$intersecting_cells, 'list'),
            expect_is(test_1$intersecting_cells[1], 'list'),
            expect_equal(length(test_1$intersecting_cells), 5),
            expect_equal(length(test_1$intersecting_cells[1]), 1),
            expect_equal(length(test_1$intersecting_cells[[1]]), 230),
            expect_is(test_1$intersecting_cells[[1]], 'integer'),
            expect_lte(max(unlist(test_1$intersecting_cells)), ncell(heronvale_covariates)),
            expect_gt(min(unlist(test_1$intersecting_cells)), 0),
            # cells don't get allocated to multiple adjacent polygons
            expect_length(base::intersect(test_1$intersecting_cells[[2]],
                                         test_1$intersecting_cells[[5]]), 0),
            # attributes carry through
            expect_identical(str(test_1[, c(1L:9L)]),
                             str(st_set_geometry(heronvale_soilmap[c(1L,2L,11L, 75L, 83L), ], NULL)))
            )
  })

test_that('prep_polygons with area_proportional rate - defaults', {
  c(
            test_2 <- prep_polygons(src_map = heronvale_soilmap[c(1L,2L,11L,75L, 83L), ],
                                    covariates    = heronvale_covariates,
                                    id_field      = 'POLY_UID',
                                    sample_method = 'area_p'),
            # default area_p sampling is default samp_floor -  2x n classes if no
            # params are provided
            # undecided whether to just leave that, its p low - only meant for polys
            # that are small relative to cell size
            expect_identical(test_2$n_soils * 2L, test_2$n_samples)
            )
  })

test_that('prep_polygons with area_proportional rate - 10/sqkm', {
  c(
    test_3 <- prep_polygons(src_map = heronvale_soilmap[c(1L,2L,11L,75L, 83L), ],
                            covariates    = heronvale_covariates,
                            id_field      = 'POLY_UID',
                            sample_method = 'area_p',
                            sample_rate = 10),
    expect_identical(test_3$n_samples, c(4L, 6L, 2L, 6L, 4L))
    )
  })

test_that('prep_polygons with area_proportional rate - no floor', {
  c(
    # force no samp_floor
    test_4 <- prep_polygons(src_map = heronvale_soilmap[c(1L,2L,11L,75L, 83L), ],
                            covariates    = heronvale_covariates,
                            id_field      = 'POLY_UID',
                            sample_method = 'area_p',
                            sample_rate = 10,
                            rate_floor = 0),
    expect_identical(test_4$n_samples, as.integer(ceiling(test_4$area_sqkm * 10)))
    )
  })

test_that('prep_polygons with area_proportional rate with cap', {
  c(
    # cap samples where only one soil class is present
    test_5 <- prep_polygons(src_map = heronvale_soilmap[c(1L,2L,11L,75L, 83L), ],
                            covariates    = heronvale_covariates,
                            id_field      = 'POLY_UID',
                            sample_method = 'area_p',
                            sample_rate = 100,
                            rate_ceiling = 5),
    expect_equal(test_5$n_samples[3], 5L)
    )
  })

#            # only n_samples should change with method
#            expect_identical(dplyr::select(test_1, -n_samples),
#                             dplyr::select(test_2, -n_samples),
#                             dplyr::select(test_3, -n_samples),
#                             dplyr::select(test_4, -n_samples),
#                             dplyr::select(test_5, -n_samples))
#          ))
#
