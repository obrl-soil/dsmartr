context('n_things')

test_that('n_things', {

  test_df <- tibble::tribble(~ID, ~C_1, ~C_2, ~P_1, ~P_2,
                               1,  'A',  'B',   90,   10,
                               2,  'A',   NA,  100,   NA,
                               3,   NA,   NA,   NA,   NA)
  expect_equal(n_things(test_df[1, ], 'C'), c('A', 'B'))
  expect_equal(n_things(test_df[1, ], 'P'), c(90, 10))
  expect_equal(n_things(test_df[2, ], 'C'), c('A'))
  expect_equal(n_things(test_df[2, ], 'P'), c(100))
  expect_true(is.na(n_things(test_df[3,], 'C')))
  expect_true(is.na(n_things(test_df[3,], 'P')))
  expect_equal(n_things(test_df, 'C'), c('A', 'A', 'B'))
  expect_equal(n_things(test_df, 'P'), c(90, 100, 10))
})
