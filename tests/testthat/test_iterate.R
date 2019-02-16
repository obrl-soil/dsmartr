context('dsmartr_iterate')

test_that('get_classes', {
  c(
    test_df  <- tibble::tribble(~ID, ~C_1, ~C_2, ~P_1, ~P_2,
                                  1,  'A',  'B',   90,   10,
                                  2,  'A',   NA,  100,   NA,
                                  3,  'D',  'C',   75,  25),
    test_pts <- tibble::tribble(~ID,   ~C,
                                  1,  'A',
                                  2,  'E',
                                  3,  'D'),
    # points data with inconsistent colname stub vs compos data
    test_pt2 <- tibble::tribble(~ID,   ~N,
                                  1,  'A',
                                  2,  'E',
                                  3,  'D'),
    expect_true(is.factor(get_classes(test_df, col_stub = 'C'))),
    expect_identical(get_classes(test_df, col_stub = 'C'),
                     factor(c('A', 'B', 'C', 'D'))),
    expect_equal(get_classes(test_df, test_pts, col_stub = 'C'),
                 factor(c('A', 'B', 'C', 'D', 'E'))),
    expect_error(get_classes(test_df, test_pt2, 'C'))
          )
  })

test_that('iter_sample_poly', {
  c(
    # from dsmartr_prep_polygons with demo dataset
    data('heronvale_soilmap'),
    data('heronvale_covariates'),
    pr_ap <- dsmartr::prep_polygons(src_map = heronvale_soilmap[c(1L,2L,11L), ],
                                    covariates = heronvale_covariates,
                                    id_field = 'POLY_UID',
                                    sample_method = 'flat',
                                    sample_rate = 100L),

    test1 <- dsmartr:::iter_sample_poly(pd = pr_ap[1L, ], cs = 'CLASS', ps = 'PERC',
                              nscol = 'n_samples', cellcol = 'intersecting_cells',
                              t_factor = 100L),
    expect_output(str(test1), '100 obs'),
    expect_output(str(test1), '2 variables'),
    expect_true(is.data.frame(test1)),
    expect_true(class(test1$CLASS) == 'character'),
    expect_false(all(is.na(match(unlist(test1$CELL),
                                 unlist(pr_ap[1, ]$intersecting_cells)))))
    )
  })


#test_that('dsmartr_iterate',
#          )





