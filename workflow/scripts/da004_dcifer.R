# DA004: Dcifer relatedness for pairs that include polyclonal samples.
# usage: Rscript da004_dcifer.R long_form.tsv coi.tsv out_pairs.tsv alpha nr lrank afreq_tol
# (layer5 da004_dcifer_alpha, da004_dcifer_null_grid, da004_dcifer_coi_lrank, da004_dcifer_afreq_tol)
suppressPackageStartupMessages({library(dcifer)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 7) stop("usage: da004_dcifer.R long_form.tsv coi.tsv out_pairs.tsv alpha nr lrank afreq_tol")
alpha <- as.numeric(a[4])
nr <- as.integer(a[5])
lrank <- as.integer(a[6])
afreq_tol <- as.numeric(a[7])
df <- read.delim(a[1], stringsAsFactors = FALSE)
coi_tab <- read.delim(a[2], stringsAsFactors = FALSE)
dsmp <- formatDat(df, svar = "sample_id", lvar = "locus", avar = "allele")
given <- setNames(suppressWarnings(as.numeric(coi_tab$coi)), coi_tab$sample_id)[names(dsmp)]
# Samples without a COI count (DA001 withheld it, or it was not run) use Dcifer's naive
# estimate from allele counts (Gerlovina et al. 2022), never an assumed 1.
naive <- getCOI(dsmp, lrank = lrank)
coi <- pmax(1L, as.integer(round(ifelse(is.na(given), naive[names(dsmp)], given))))
names(coi) <- names(dsmp)
cat("COI from DA001:", sum(!is.na(given)), "; naive getCOI:", sum(is.na(given)), "\n")
afreq <- calcAfreq(dsmp, coi, tol = afreq_tol)
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
