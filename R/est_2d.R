#' Estimate 1 diffusion coefficient
#'
#' @description Provides the negative log likelihood used for estimating diffusion coefficient
#'
#' @param log_D1 the log of the diffusion coefficient for movement state 1
#' @param log_D2 the log of the diffusion coefficient for movement state 2
#' @param h Resolution, or the length of a side in meters
#' @param L Likelihood for your fish over space and indexed by time
#' @param land Raster of land so that any probability on land can be set to zero each step
#' @param fish_data Data on individual fish, used here for multiple movement states (mvst)
#'
#' @return numeric value of the negative log likelihood
#' @export
#'
#' @examples est_1d(log_D = log(70), h = h, L = L, land = land)
est_2d <- function(log_D1 = log(D1), log_D2 = log(D2), h = h, L = L, land = land, fish_data = fish_data){
  D1 <- exp(log_D1)
  D2 <- exp(log_D2)
  sig.kern1 <- sqrt(2 * D1 / ((h / 1000) ^ 2)) # h/1000 converts D into map units, in km
  kern1 <- conv_kern(sigma = sig.kern1, k = 'gaussian')
  sig.kern2 <- sqrt(2 * D2 / ((h / 1000) ^ 2)) # h/1000 converts D into map units, in km
  kern2 <- conv_kern(sigma = sig.kern2, k = 'gaussian')
  pred <- array(0, dim = dim(L)) #predicted
  post <- L[, , 1] #Location at Day one is posterior for Day 1
  icalc <- dim(L)[3] # number of days at liberty
  lambda <- rep(NA, icalc - 1) #Holds sum of probabilities for each time step, basis of maximu likelihood

  for (i in 2:icalc) {
    if (fish_data$mvst[i] == 1) {
      kern <- kern1
    }
    else if (fish_data$mvst[i] == 2) {
      kern <- kern2
    }
    p1 <- imager::as.cimg(t(post))
    K <- imager::as.cimg(kern$matrix)
    q <- imager::convolve(p1, K)
    q1 <- t(terra::as.matrix(q))
    q1[land == 0] <- 0
    pred[, , i] <- q1
    post <- L[, , i] * q1
    lambda[i - 1] <- sum(post)
    post <- post / sum(post)
  }
  nll <- -1 * sum(log(lambda)) #Negative log likelihood
  return(nll)
}
