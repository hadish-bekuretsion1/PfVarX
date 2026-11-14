# DA004: Dcifer relatedness for pairs that include polyclonal samples.
# usage: Rscript da004_dcifer.R long_form.tsv coi.tsv out_pairs.tsv [alpha] [nr]
suppressPackageStartupMessages({library(dcifer)})
a <- commandArgs(trailingOnly = TRUE)
alpha <- if (length(a) >= 4) as.numeric(a[4]) else 0.05
nr <- if (length(a) >= 5) as.integer(a[5]) else 1000L
df <- read.delim(a[1], stringsAsFactors = FALSE)
coi_tab <- read.delim(a[2], stringsAsFactors = FALSE)
dsmp <- formatDat(df, svar = "sample_id", lvar = "locus", avar = "allele")
coi <- setNames(pmax(1L, as.integer(round(coi_tab$coi))), coi_tab$sample_id)[names(dsmp)]
afreq <- calcAfreq(dsmp, coi, tol = 1e-5)
res <- ibdDat(dsmp, coi, afreq, pval = TRUE, confint = TRUE, rnull = 0, alpha = alpha, nr = nr)
ids <- names(dsmp)
out <- data.frame()
for (i in seq_along(ids)) for (j in seq_along(ids)) if (i < j) {
  out <- rbind(out, data.frame(sample_1 = ids[j], sample_2 = ids[i],
    estimate = res[j, i, "estimate"], ci_lower = res[j, i, "CI_lower"],
    ci_upper = res[j, i, "CI_upper"], p_value = res[j, i, "p_value"]))
}
write.table(out, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("dcifer", as.character(packageVersion("dcifer")), "\n")
