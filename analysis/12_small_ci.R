#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 12_small_ci.R  (any machine, no network)
#
# Runs the 800-cell subsample twice, under both interval estimators, and writes
# the numbers the Supplementary quotes for the subsample: the medians, the
# delete-one-donor jackknife interval, and the bootstrap interval that the
# jackknife replaced. Kept separate from 11_small_reproduce.R so that the
# shipped facts file (small_dataset.txt) stays free of the comparison with the
# superseded estimator.
#
# Usage:
#   Rscript 12_small_ci.R [pkgdir] [outdir]
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args) >= 1) args[1] else "."
outdir <- if (length(args) >= 2) args[2] else tempdir()

suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

rds <- file.path(pkg, "inst", "extdata", "gse96583_small.rds")
stopifnot(file.exists(rds))
small <- readRDS(rds)
stopifnot(all(c("counts", "donor", "group") %in% names(small)))

tb <- table(small$donor, small$group)
stopifnot(ncol(tb) == 2L, all(tb == 50L))

run <- function(jack) {
  effn(small$counts, group = small$group, donor = small$donor,
       n_boot = 200L, n_tech = 40L, seed = 20261021L, jack = jack)
}

jk <- run(TRUE)
bs <- run(FALSE)

stopifnot(identical(jk$params$jack, TRUE), identical(bs$params$jack, FALSE))
stopifnot(identical(jk$summary$ci_method, "jackknife (delete one donor)"))

sm <- function(x) x$summary
lines <- c(
  sprintf("median_mde_log2=%.10f", sm(jk)$median_mde_log2),
  sprintf("median_n_eff=%.10f", sm(jk)$median_n_eff),
  sprintf("jack_lo=%.10f", sm(jk)$median_mde_log2_ci[1]),
  sprintf("jack_hi=%.10f", sm(jk)$median_mde_log2_ci[2]),
  sprintf("boot_lo=%.10f", sm(bs)$median_mde_log2_ci[1]),
  sprintf("boot_hi=%.10f", sm(bs)$median_mde_log2_ci[2]),
  sprintf("jack_contains_point=%s",
          sm(jk)$median_mde_log2_ci[1] <= sm(jk)$median_mde_log2 &&
            sm(jk)$median_mde_log2 <= sm(jk)$median_mde_log2_ci[2]),
  sprintf("boot_contains_point=%s",
          sm(bs)$median_mde_log2_ci[1] <= sm(bs)$median_mde_log2 &&
            sm(bs)$median_mde_log2 <= sm(bs)$median_mde_log2_ci[2]),
  sprintf("n_donor=%d", length(unique(small$donor)))
)
writeLines(lines, file.path(outdir, "small_dataset_ci.txt"))
cat(lines, sep = "\n")
cat("\nwrote", file.path(outdir, "small_dataset_ci.txt"), "\n")
cat("SMALLCIDONE\n")
