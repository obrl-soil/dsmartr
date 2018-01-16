context('collate')

test_that('count_predictions', {
  c(
    cell_19291 <- as.integer(
      c(29, 7, 7, 7, 12, 9, 29, 29, 29, 3, 16, 15, 9, 21, 24, 2, 9, 15, 34, 9, 24, 34, 24, 34, 21,
        24, 15, 15, 9, 21, 21, 7, 7, 34, 34, 15, 15, 21, 7, 30, 32, 21, 29, 21, 24, 29, 31, 12, 24,
        12, 21, 15, 15, 29, 12, 7, 21, 12, 7, 34, 9, 32, 7, 7, 32, 32, 12, 34, 29, 7, 3, 15, 15, 29,
        7, 3, 34, 9, 9, 33, 32, 21, 34, 29, 33, 7, 7, 31, 7, 7, 24, 21, 15, 7, 24, 29, 21, 31, 29,
        15)),
    cell_null <- rep(NA_integer_, 100),
    cell_mix  <- cell_19291,
    # ponder this some more
    cell_mix[c(5,25,50,75)] <- NA_integer_,
    expect_identical(count_predictions(cell_19291, 38),
                     as.integer(c( 0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0,
                                  12, 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0))),
    expect_identical(count_predictions(cell_null, 38), rep(NA_integer_, 38)),
    expect_identical(count_predictions(cell_mix, 38),
                     as.integer(c( 0, 1, 3, 0,  0, 0, 16,  0, 8, 0, 0, 4, 0, 0, 12, 1, 0, 0, 0, 0,
                                  11, 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)))
  )
})

test_that('calc_probabilities', {
  c(
    counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
                      0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    expect_identical(calc_probabilities(counts_19291, 100)[7], 0.170),
    expect_identical(calc_probabilities(rep(NA_integer_, 38), 100), rep(NA_real_, 38))
  )
})

test_that('order_counts', {
  c(
    counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
                                 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    expect_identical(order_counts(counts_19291, 38)[1], 7L),
    expect_identical(order_counts(counts_19291, 38)[5], 34L),
    # shuffle test
    expect_setequal(order_counts(counts_19291, 38)[2:4], c(15L, 21L, 29L))
  )
})

test_that('sort_probabilities', {
  c(
    counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
                                 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    probs_19291 <- calc_probabilities(counts_19291, 100),
    expect_identical(sort_probabilities(probs_19291, 100),
                     c(0.170, 0.120, 0.120, 0.120, 0.090, 0.080, 0.080, 0.060, 0.050, 0.030, 0.030,
                       0.020, 0.010, 0.010, 0.010, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                       0.000, 0.000, 0.000, 0.000, 0.000))
  )
})

