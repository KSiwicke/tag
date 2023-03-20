est_1d <- function(log_D = log(D)){
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
    q1 <- t(as.matrix(q))
    q1[land==0] <- 0
    pred[, , i] <- q1
    post <- L[, , i] * q1
    lambda[i - 1] <- sum(post)
    post <- post / sum(post)
  }
  nll <- -1 * sum(log(lambda)) #Negative log likelihood
  return(nll)
}
