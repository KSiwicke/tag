#' Reverse Latitudes
#' @description Reverse latitudes for plotting matrix as an image
#'
#' @param m a matrix of latitude values
#'
#' @return a matrix of reversed latitude values
#' @export
#'
#' @examples
#' m <- matrix(c(50, 55, 40, 45), nrow = 2)
#' rev_lat(m)
rev_lat <- function(m) {
  t(m)[, nrow(m):1]
}
