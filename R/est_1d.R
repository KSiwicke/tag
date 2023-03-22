#' Estimate 1 diffusion coefficient
#'
#' @description Provides the negative log likelihood used for estimating diffusion coefficient
#'
#' @param log_D the log of the diffusion coefficient
#' @param h Resolution, or the length of a side in meters
#' @param L Likelihood for your fish over space and indexed by time
#' @param land Raster of land so that any probability on land can be set to zero each step
#'
#'
#' @return numeric value of the negative log likelihood
#' @export
#'
#' @examples est_1d(log_D = log(70), h = h, L = L, land = land)
est_1d <- function(log_D = log(D), h = h, L = L, land = land){
  D <- exp(log_D)
  sig <- sqrt(2 * D / (h / 1000)^2) # resolution /1000 converts D into map units, in km
  kern <- conv_kern(sigma = sig, k = 'gaussian')
  pred <- array(0, dim = dim(L)) #predicted
  post <- L[, , 1] #Location at Day one is posterior for Day 1
  lambda <- rep(NA, length(t) - 1) #Holds sum of probabilities for each time step, basis of maximu likelihood

  for (i in 2:length(t)) {
    p1 <- imager::as.cimg(t(post))
    K <- imager::as.cimg(kern$matrix)
    q <- imager::convolve(p1, K)
    q1 <- t(terra::as.matrix(q))
    q1[land==0] <- 0
    pred[, , i] <- q1
    post <- L[, , i] * q1
    lambda[i - 1] <- sum(post)
    post <- post / sum(post)
  }
  nll <- -1 * sum(log(lambda)) #Negative log likelihood
  return(nll)
}
