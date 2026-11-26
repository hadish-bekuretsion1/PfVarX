# DA007: binomial GAM of marker prevalence on time, smoothness by REML (mgcv).
# usage: Rscript da007_gam.R input.tsv out.tsv max_basis interval_level grid_points
# (layer5 da007_gam_max_basis, da007_gam_interval_level, da007_gam_grid_points)
suppressPackageStartupMessages({library(mgcv)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 5) stop("usage: da007_gam.R input.tsv out.tsv max_basis interval_level grid_points")
kmax <- as.integer(a[3])
z <- qnorm(1 - (1 - as.numeric(a[4])) / 2)
npoints <- as.integer(a[5])
d <- read.delim(a[1])
out <- data.frame()
for (m in unique(d$marker_id)) {
  x <- d[d$marker_id == m, ]
  if (length(unique(x$t)) < kmax || length(unique(x$y)) < 2) next
  k <- min(kmax, length(unique(x$t)) - 1)
  fit <- gam(y ~ s(t, k = k), family = binomial, data = x, method = "REML")
  grid <- data.frame(t = seq(min(x$t), max(x$t), length.out = npoints))
  p <- predict(fit, grid, se.fit = TRUE)
  out <- rbind(out, data.frame(marker_id = m, t = grid$t, fit = plogis(p$fit),
                               lo = plogis(p$fit - z * p$se.fit), hi = plogis(p$fit + z * p$se.fit),
                               edf = sum(fit$edf)))
}
write.table(out, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
