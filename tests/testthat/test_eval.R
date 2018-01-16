context('eval')

test_that('n_predicted', {
  c(
    counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
                                 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    expect_identical(n_predicted(counts_19291), 15L),
    expect_identical(n_predicted(counts_19291, 100, 0.1), 4L),
    expect_error(n_predicted(counts_19291, 100, 10))
  )
})

test_that('tie_finder', {
  c(
    counts_19291 <- as.integer(c(0, 1, 3, 0,  0, 0, 17,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0, 12,
                                 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    counts_19291_mod <- as.integer(c(0, 1, 3, 0,  0, 0, 6,  0, 8, 0, 0, 6, 0, 0, 12, 1, 0, 0, 0, 0,
                                     12, 0, 0, 8, 0,  0, 0,  0, 12, 1, 3, 5, 2, 9, 0,  0, 0, 0)),
    expect_identical(tie_finder(counts_19291), 0L),
    expect_identical(tie_finder(counts_19291_mod), 3L),
    expect_identical(tie_finder(rep(NA_integer_, 38)), NA_integer_)
  )
})

