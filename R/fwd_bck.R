#' Forward filter for 1- or 2-state movement to estimate diffusion.
#'
#' @param D Diffusion coefficient for 1 or 2 movement states, can be preselected or estimated
#' @param h Resolution, or the length of a side in meters
#' @param L Likelihood for your fish over space and indexed by time
#' @param fish.data Data on individual fish, used here for multiple movement states (mvst)
#' @param land Raster of land so that any probability on land can be set to zero each step
#'
#' @return A numeric negaitve-log-likelihood value.
#' @export
#'
#' @examples
#' fwd_filter(50.3)
#' est.D <- optim(50, fwd_filter, lower = 2, upper = 300, method = "Brent")
#' # need to choose upper and lower boundaries

fwd_bck <- function(D = D, h = h, L = L, fish.data = fish.data, land = land) {
  pred <- array(0, dim = dim(L)) # predicted
  phi <- array(0, dim = dim(L))  # holds probability
  phi[, , 1] <- L[, , 1] # All probability in release location, dirac delta
  post <- L[, , 1] # Location at Day one is posterior for Day 1
  icalc <- dim(L)[3] # number of days at liberty
  # lambda <- rep(NA, icalc - 1) # Holds sum of probabilities for each time step, basis of maximum likelihood
  sig.kern1 <- sqrt(2 * D[1] / ((h / 1000) ^ 2)) # h/1000 converts D into map units, in km
  kern1 <- conv_kern(sigma = sig.kern1, k = 'gaussian')

  if (length(D) > 1) {
    sig.kern2 <- sqrt(2 * D[2] / ((h / 1000) ^ 2)) # h/1000 converts D into map units, in km
    kern2 <- conv_kern(sigma = sig.kern2, k = 'gaussian')
  }

  for (i in 2:icalc) {
    if (fish.data$mvst[i] == 1) {
      kern <- kern1
    }
    else if (fish.data$mvst[i] == 2) {
      kern <- kern2
    }
    p1 <- imager::as.cimg(t(post))
    K <- imager::as.cimg(kern$matrix)
    q <- imager::convolve(p1, K)
    q1 <- t(terra::as.matrix(q))
    q1[land == 0] <- 0
    pred[, , i] <- q1
    post <- L[, , i] * q1
    # lambda[i - 1] <- sum(post)
    post <- post / sum(post)
    phi[, , i] <- post
  }
  # Backwards smoother
  smooth <- array(0, dim = dim(L))
  smooth[, , icalc] <- phi[, , icalc] #Start at the end and work toward beginning

  for(i in icalc:2) {
    ratio <- smooth[, , i] / (pred[, , i] + 1e-15)
    p1 = imager::as.cimg(t(ratio))
    if (fish.data$mvst[i] == 1) {
      kern <- kern1
    } else if (fish.data$mvst[i] == 2) {
      kern <- kern2
    }
    K <- imager::as.cimg(kern$matrix)
    Rp1 <- imager::convolve(p1, K)
    Rp1 = t(terra::as.matrix(Rp1))
    Rp1[land == 0] <- 0 # Added that was not present before...KAS
    smooth[, , i - 1] <- phi[, , i - 1] * Rp1
    smooth[, , i - 1] <- smooth[, , i - 1] / sum(smooth[, , i - 1])
  }
  return(smooth)
}
