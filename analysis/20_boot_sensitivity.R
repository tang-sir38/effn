#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 20_boot_sensitivity.R
# How many technical resamples does the depth term need?
#
# n_tech is the number of times each donor's cells are resampled to estimate
# sigma_tech^2, so it is the knob that sets the Monte Carlo error of the whole
# decomposition. The default is 40. This script runs the shipped 800-cell
# subsample at five values of n_tech and three seeds each, so that both the
# bias (distance from the largest n_tech) and the noise (spread across seeds)
# are visible. n_boot is left at 200: it drives the confidence interval, not
# the point estimate, and is not the quantity under test here.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(effn))

SM <- Sys.getenv("EFFN_SMALL", "C:/rtmp/effnfig2/gse96583_small.rds")
OUT <- Sys.getenv("EFFN_FIG", "C:/rtmp/effnfig2")
s <- readRDS(SM)
cat("subsample:", nrow(s$counts), "genes x", ncol(s$counts), "cells\n")

NTECH <- c(5L, 10L, 20L, 40L, 80L)
SEEDS <- c(20261021L, 20261022L, 20261023L)

rows <- list()
for (b in NTECH) {
  for (sd in SEEDS) {
    t0 <- Sys.time()
    r <- effn(s$counts, group = s$group, donor = s$donor,
              n_boot = 200L, n_tech = b, seed = sd)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok <- r$per_gene[r$per_gene$status == "ok", ]
    rows[[length(rows) + 1L]] <- data.frame(
      n_tech = b, seed = sd,
      median_sigma_tech = stats::median(ok$sigma_tech, na.rm = TRUE),
      median_n_eff = r$summary$median_n_eff,
      median_mde_fold = r$summary$median_mde_fold,
      median_rho = r$summary$median_rho_tech,
      n_ok = nrow(ok), secs = el,
      stringsAsFactors = FALSE)
    cat(sprintf("n_tech=%3d seed=%d  sigma_tech=%.6f  n_eff=%.4f  mde=%.4f  %.1fs\n",
                b, sd, utils::tail(rows, 1)[[1]]$median_sigma_tech,
                r$summary$median_n_eff, r$summary$median_mde_fold, el))
  }
}
res <- do.call(rbind, rows)
utils::write.table(res, file.path(OUT, "t_boot_sensitivity.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

## summary per n_tech: bias against the largest n_tech, noise across seeds
ref <- mean(res$median_sigma_tech[res$n_tech == max(NTECH)])
agg <- do.call(rbind, lapply(NTECH, function(b) {
  d <- res[res$n_tech == b, ]
  data.frame(n_tech = b,
             sigma_tech_mean = mean(d$median_sigma_tech),
             bias_pct = 100 * (mean(d$median_sigma_tech) - ref) / ref,
             seed_cv_pct = 100 * stats::sd(d$median_sigma_tech) /
               mean(d$median_sigma_tech),
             mde_mean = mean(d$median_mde_fold),
             secs_mean = mean(d$secs))
}))
utils::write.table(agg, file.path(OUT, "t_boot_sensitivity_summary.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(agg, digits = 3)
cat("BOOTDONE\n")
