Dsim = 70
h = 3000
days = 60
bathy <- cod_ex[[1]]
st_dimensions(bathy)$x$delta # resolution in x direction
st_dimensions(bathy)$y$delta # resolution in y direction
h <- st_dimensions(bathy)$x$delta

land <- bathy %>%
  mutate(depth = ifelse(depth < 15, 0, 1)) # 15-m cutoff, but can change

land = t(land$depth)

start = c(round(runif(1, min = 0, max = nrow(bathy))), round(runif(1, min = 0, max = ncol(bathy))))

sim_Loc <- function(D = Dsim, bathy = bathy, h = h, start = start, days = days) {
  sig <- sqrt(2 * D / (h / 1000)^2) # resolution /1000 converts D into map units, in km
  kern <- conv_kern(sigma = sig, k = 'gaussian')
  icalc <- days # number of days at liberty
  L <- array(0, dim = c(dim(bathy), days))
  L[start[1], start[2], 1] <- 1
  index <- matrix(1:length(land), nrow = nrow(land))

  for (i in 2:icalc) {
    p1 <- imager::as.cimg(t(L[, , i-1]))
    K <- imager::as.cimg(kern$matrix)
    q <- imager::convolve(p1, K)
    q1 <- t(terra::as.matrix(q))
    q1[land == 0] <- 0
    probs <- q1 / sum(q1)

    loc <- sample(size = 1, x = 1:length(probs), prob = as.vector(probs))
    L[, , i] <- ifelse( index == loc, 1, 0)
  }
}
