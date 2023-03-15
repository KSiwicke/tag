test_that("convolution kernel works", {
  expect_equal(conv_kern(sigma = 0.5, k = 'gaussian')$matrix,
               matrix(c(0.011660098, 0.086157117, 0.011660098, 0.086157117, 0.636619772, 0.086157117, 0.011660098, 0.086157117, 0.011660098), nrow = 3)
  )
})
