test_that("rev_lat reverses latitude values", {
  expect_equal(rev_lat(matrix(c(50, 55, 40, 45), nrow = 2)), matrix(c(55, 45, 50, 40), nrow = 2))
})
