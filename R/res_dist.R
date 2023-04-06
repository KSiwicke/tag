res_dist <- function(smooth = smooth) {
  smooth.summary <- apply(smooth, c(1, 2), sum) #Sum over entire time period
  smooth.summary.norm <- smooth.summary / sum(smooth.summary)
  res.distr <- cbind(1:length(smooth.summary.norm), as.vector(smooth.summary.norm))
  res.distr.order <- res.distr[order(smooth.summary.norm, decreasing = F), ]
  res.distr.cumsum <- cbind(res.distr.order, cumsum(res.distr.order[, 2]))
  res.distr.quantile <- matrix(NA, rows,cols)
  res.distr.quantile[res.distr.cumsum[,1]] <- res.distr.cumsum[, 3]
  image.plot(longs + h / 2, lats + h / 2, reverse(res.distr.quantile), col = seq.pal(100), main="Residence distribution", xlab = "", ylab = "")
  image(longs + h / 2, lats + h / 2, reverse(land), zlim = c(0, 0.001), add = T, col = "black")
}

