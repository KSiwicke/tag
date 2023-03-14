rev_lat <- function(m) {
  t(m)[,nrow(m):1]
}
