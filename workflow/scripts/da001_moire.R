# DA001: MOIRE (Murphy & Greenhouse) on long-form allele-presence data.
# usage: Rscript da001_moire.R long_form.tsv out_coi.tsv out_afreq.tsv burnin samples seed
suppressPackageStartupMessages({library(moire)})
a <- commandArgs(trailingOnly = TRUE)
df <- read.delim(a[1], stringsAsFactors = FALSE)
dat <- load_long_form_data(df, warn_uninformative = FALSE)
res <- run_mcmc(dat, is_missing = dat$is_missing, burnin = as.integer(a[4]),
                samples_per_chain = as.integer(a[5]), verbose = FALSE, seed = as.integer(a[6]),
                pt_chains = 1)
coi <- summarize_coi(res)
af <- summarize_allele_freqs(res)
write.table(coi, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(af, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("moire", as.character(packageVersion("moire")), "\n")
