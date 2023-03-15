#' Create a convolution kernel
#'
#' @description This function creates the convolution kernel for applying a filter to an array/matrix
#' - spatialfil package is not available anymore, need to load function manually
#'
#' @param sigma The numeric value of standard deviation for the Gaussian or LoG filter
#' @param k character value:
#' "gaussian" for Gaussian kernel,
#' "LoG" for Laplacian of Gaussian kernel,
#' "sharpen" for 3x3 convolution matrix for sharpening edges, or
#' "laplacian" for a 3x3 convolution matrix.
#'
#' @return A numeric convolution matrix
#' @export
#'
#' @examples # Not run:
#' # creates a convolution kernel with Gaussian function and sigma = 1.4
#' K <- conv_kern(sigma = 1.4, k = 'gaussian')
conv_kern <- function(sigma = 1.4, k = c("gaussian", "LoG", "sharpen",
                                          "laplacian", "emboss", "sobel")) {
  k <- match.arg(k)
  l <- sigma * 7
  if (l%%2 == 0)
    l <- l + 1
  x <- c(-floor(l / 2):floor(l / 2))
  y <- c(-floor(l / 2):floor(l / 2))
  if (k == "gaussian")
    M <- outer(X = x, Y = y, FUN = function(X, Y)
      return(1 / (2 * pi * sigma ^ 2) * exp(-(X ^ 2 + Y ^ 2) / (2 * sigma ^ 2))))
  if (k == "LoG")
    M <- outer(X = x, Y = y, FUN = function(X, Y)
      return(-1 / (pi * sigma ^ 4) * (1 - (X ^ 2 + Y ^ 2) / (2 * sigma ^ 2)) *
               exp(-(X ^ 2 + Y ^ 2) / (2 * sigma ^ 2))))
  if (k == "sharpen")
    M <- matrix(data = c(0, -1, 0, -1, 5, -1, 0, -1, 0), nrow = 3)
  if (k == "laplacian")
    M <- matrix(data = c(0.5, 1, 0.5, 1, -6, 1, 0.5, 1, 0.5), nrow = 3)
  if (k == "emboss")
    M <- matrix(c(2, 0, 0, 0, -1, 0, 0, 0, -1), nrow = 3)
  if (k == "sobel")
    M <- matrix(c(1, 2, 1, 0, 0, 0, -1, -2, -1), nrow = 3)
  if ((k == "LoG") || (k == "gaussian"))
    output <- list(matrix = M, kernel = k, sigma = sigma)
  else output <- list(matrix = M, kernel = k)
  class(output) <- "convKern"
  return(output)
}
