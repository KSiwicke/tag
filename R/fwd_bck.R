#' Forward filter for 1- or 2-state movement to estimate diffusion.
#'
#' @param D Diffusion coefficient for 1 or 2 movement states, can be preselected or estimated
#' @param h Resolution, or the length of a side in meters
#' @param L Likelihood for your fish over space and indexed by time
#' @param fish_data Data on individual fish, used here for multiple movement states (mvst)
#' @param land Raster of land so that any probability on land can be set to zero each step
#'
#' @return A numeric negaitve-log-likelihood value.
#' @export
#'
#' @examples
#' h <- 1000 # 1 km resolution
#' set.seed(43516)
#' bathy <- matrix(rnorm(100, 50, 30), nrow = 10)
#' land <- ifelse(bathy > 0, 1, 0)
#' L <- array(0, dim = c(dim(bathy), 2))
#' L[5, 5, 1] <- 1
#' L[6, 6, 2] <- 1
#' fish_data <- data.frame("mvmt" = c(1, 1))
#' fwd_bck(D = 10, h = h, L = L, land = land, fish_data = fish_data)

fwd_bck <- function(D = D, h = h, L = L, fish_data = fish_data, land = land) {
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
    kern <- kern1
    if ( length(D) > 1) {
      if (fish_data$mvst[i] == 2) {
        kern <- kern2
      }
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
    kern <- kern1
    if ( length(D) > 1) {
      if (fish_data$mvst[i] == 2) {
        kern <- kern2
      }
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
