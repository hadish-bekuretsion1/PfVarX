# DA007: binomial GAM of marker prevalence on time, smoothness by REML (mgcv).
# usage: Rscript da007_gam.R input.tsv out.tsv
suppressPackageStartupMessages({library(mgcv)})
a <- commandArgs(trailingOnly = TRUE)
d <- read.delim(a[1])
out <- data.frame()
for (m in unique(d$marker_id)) {
  x <- d[d$marker_id == m, ]
  if (length(unique(x$t)) < 4 || length(unique(x$y)) < 2) next
  k <- min(4, length(unique(x$t)) - 1)
  fit <- gam(y ~ s(t, k = k), family = binomial, data = x, method = "REML")
  grid <- data.frame(t = seq(min(x$t), max(x$t), length.out = 20))
  p <- predict(fit, grid, se.fit = TRUE)
  out <- rbind(out, data.frame(marker_id = m, t = grid$t, fit = plogis(p$fit),
                               lo = plogis(p$fit - 1.96 * p$se.fit), hi = plogis(p$fit + 1.96 * p$se.fit),
                               edf = sum(fit$edf)))
}
write.table(out, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
