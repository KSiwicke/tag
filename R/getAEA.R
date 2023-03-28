#' Convert from latitude/longitude to Albers Equal Area, same as Aleutians map (bathy.ai)

#'
#' @description Provides the negative log likelihood used for estimating diffusion coefficient
#'
#' @param lon longitude
#' @param lat latitude
#'
#' @return x and y coordinates in Albers Equal Area with names 'X' and 'Y
#' @export
getAEA <- function(lon, lat) {
  xy <- terra::project(cbind(lon, lat), proj = "+proj=aea +lat_0=50 +lon_0=-154 +lat_1=55 +lat_2=65 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")
  colnames(xy) <- c("X", "Y")
  return(xy)
}
