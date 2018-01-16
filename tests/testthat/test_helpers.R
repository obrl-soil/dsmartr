context('helpers')

test_that('n_things', {
  c(
  test_df <- tibble::tribble(~ID,            ~C_1,           ~C_2,          ~P_1,          ~P_2,
                               1,             'A',            'B',           90L,           10L,
                               2,             'A',  NA_character_,          100L,   NA_integer_,
                               3,   NA_character_,  NA_character_,   NA_integer_,   NA_integer_ ),
  expect_equal(n_things(test_df[1, ], 'C'), c('A', 'B')),
  expect_equal(n_things(test_df[1, ], 'P'), c(90L, 10L)),
  expect_equal(n_things(test_df[2, ], 'C'), c('A')),
  expect_equal(n_things(test_df[2, ], 'P'), c(100L)),
  expect_true(is.na(n_things(test_df[3,], 'C'))),
  expect_true(is.na(n_things(test_df[3,], 'P'))),
  expect_equal(n_things(test_df, 'C'), c('A', 'A', 'B')),
  expect_equal(n_things(test_df, 'P'), c(90L, 100L, 10L))
  )
})

test_that('in_range', {
  c(
    expect_true(in_range(1,0,2)),
    expect_true(in_range(1,-1,2)),
    expect_true(in_range(-1,-1,0, strict = FALSE)),
    expect_true(in_range(0,0,1, strict = FALSE)),
    expect_true(in_range(5,3,5, strict = FALSE)),
    expect_true(in_range(0.0001, 0, 1)),
    expect_true(in_range(0.001, 0.00001, 0.01)),
    expect_error(in_range(0,1,0))
  )
})

test_that('check_attributes', {
  c(
    data('heronvale_soilmap'),
    test_1 <- check_attributes(src_map = heronvale_soilmap, id_field = 'POLY_UID',
                               cs = 'CLASS', ps = 'PERC'),

    expect_is(test_1, 'sf'),
    expect_identical(ncol(heronvale_soilmap) + 4L, ncol(test_1)),
    expect_identical(heronvale_soilmap$geometry, test_1$geometry),
    expect_is(test_1$missing_data, 'logical'),
    expect_is(test_1$zero_percs, 'logical'),
    expect_is(test_1$problem_percs, 'logical'),
    expect_is(test_1$duplicate_ids, 'logical')
    )
  })

test_that('strict_cfp_total_intersect', {
  c(
    tpoly <- st_sf('ID'   = 1,
                   'geom' = st_sfc(st_polygon(x = list(matrix(c(2,2,2,5,5,5,5,2,2,2),
                                                              ncol = 2, byrow = TRUE))))),
    trast <- raster::raster(matrix(sample(c(1:49), size = 49, replace = FALSE), ncol = 7, byrow = T),
                            xmn = 0, xmx = 7, ymn = 0, ymx = 7, crs = NA, template = NULL),
    expect_is(strict_cfp(tpoly, trast), 'list'),
    expect_is(strict_cfp(tpoly, trast)[[1]], 'integer'),
    expect_identical(strict_cfp(tpoly, trast)[[1]],
                     as.integer(c(17, 18, 19, 24, 25, 26, 31, 32, 33)))
    )
})


test_that('strict_cfp_partial_intersect', {
  c(
    tpoly <- st_sf('ID'   = 1,
                   'geom' = st_sfc(st_polygon(x = list(matrix(c(2,2,2,5,5,5,5,2,2,2),
                                                              ncol = 2, byrow = TRUE))))),
    trast <- raster::raster(matrix(sample(c(1:32), size = 32, replace = FALSE), ncol = 4, byrow = T),
                            xmn = 0, xmx = 4, ymn = 0, ymx = 8, crs = NA, template = NULL),
    expect_is(strict_cfp(tpoly, trast), 'list'),
    expect_is(strict_cfp(tpoly, trast)[[1]], 'integer'),
    expect_identical(strict_cfp(tpoly, trast)[[1]],
                     as.integer(c(15,16,19,20,23,24)))
  )
})

test_that('strict_cfp_no_intersect', {
  c(
    tpoly <- st_sf('ID'   = 1,
                   'geom' = st_sfc(st_polygon(x = list(matrix(c(2,2,2,5,5,5,5,2,2,2),
                                                              ncol = 2, byrow = TRUE))))),
    trast <- raster::raster(matrix(sample(c(1:25), size = 25, replace = FALSE), ncol = 5, byrow = T),
                            xmn = 6, xmx = 11, ymn = 6, ymx = 11, crs = NA, template = NULL),
    expect_is(strict_cfp(tpoly, trast), 'list'),
    expect_is(strict_cfp(tpoly, trast)[[1]], 'integer'),
    expect_identical(strict_cfp(tpoly, trast)[[1]], NA_integer_)
    )
})
