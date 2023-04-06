#' Make error grid for smoothed probability of location
#'
#' @param smooth smoothed daily output of forward backward filter
#' @param fish_data Data on individual fish, used here for multiple movement states (mvst)
#'
#' @return A list for daily error grids
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
err_grid <- function(smooth = smooth, fish_data = fish_data) {
  `%do%` <- foreach::`%do%`
  err.grid <- foreach::foreach(m = 1:dim(smooth)[3], .combine = 'c', .packages = c("terra", "lubridate")) %do% {
    tmp <- smooth[, , m]
    tmp.norm <- tmp / sum(tmp) # sum of tmp is always 1, so tmp = tmp.norm
    tmp.distr <- cbind(1:length(tmp.norm), as.vector(tmp.norm)) # length same for all
    tmp.distr.order <- tmp.distr[order(tmp.norm, decreasing = F), ]
    tmp.distr.cumsum <- cbind(tmp.distr.order, cumsum(tmp.distr.order[, 2]))
    tmp.distr.quantile <- matrix(NA, dim(smooth)[1], dim(smooth)[2])
    tmp.distr.quantile[tmp.distr.cumsum[, 1]] <- tmp.distr.cumsum[, 3]
    tmp.distr.quantile.rst <- terra::rast(tmp.distr.quantile, extent = terra::ext(bathy), crs = terra::crs(bathy))
    err.grid.smooth.99 <- as.polygons(tmp.distr.quantile.rst >= .01, dissolve = TRUE)
    err.grid.smooth.99[["Date"]] <- fish_data$day[m]
    err.grid.smooth.99[["month"]] <- lubridate::month(fish_data$day[m])
    err.grid.smooth.99[["quantile"]] <- 99
    err.grid.smooth.50 <- as.polygons(tmp.distr.quantile.rst >= .5, dissolve = TRUE)
    err.grid.smooth.50[["Date"]] <- fish_data$day[m]
    err.grid.smooth.50[["month"]] <- lubridate::month(fish_data$day[m])
    err.grid.smooth.50[["quantile"]] <- 50
    output <- list(q99 = err.grid.smooth.99, q50 = err.grid.smooth.50)
  }
  return(err.grid)
}
