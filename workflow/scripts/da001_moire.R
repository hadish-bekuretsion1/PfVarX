# DA001: MOIRE (Murphy & Greenhouse) on long-form allele-presence data.
# usage: Rscript da001_moire.R long_form.tsv out_coi.tsv out_afreq.tsv burnin samples seed pt_chains
# (layer5 da001_moire_burnin, da001_moire_samples, config random_seed, da001_moire_pt_chains)
suppressPackageStartupMessages({library(moire)})
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 7) stop("usage: da001_moire.R long_form.tsv out_coi.tsv out_afreq.tsv burnin samples seed pt_chains")
df <- read.delim(a[1], stringsAsFactors = FALSE)
dat <- load_long_form_data(df, warn_uninformative = FALSE)
res <- run_mcmc(dat, is_missing = dat$is_missing, burnin = as.integer(a[4]),
                samples_per_chain = as.integer(a[5]), verbose = FALSE, seed = as.integer(a[6]),
                pt_chains = as.integer(a[7]))
coi <- summarize_coi(res)
af <- summarize_allele_freqs(res)
write.table(coi, a[2], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(af, a[3], sep = "\t", quote = FALSE, row.names = FALSE)
cat("moire", as.character(packageVersion("moire")), "\n")
