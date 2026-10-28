#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 21_nboot_sensitivity.R
# Does the bootstrap budget move the answer?
#
# n_boot is the number of donor resamples used for the confidence interval.
# The manual claim is that it does not touch the point estimate, because the
# interval and the estimate come from different loops. This script tests that
# claim rather than asserting it: five values of n_boot (10 to 400), three
# seeds each, at fixed n_tech = 40.
#
# The interval is only built from n_boot when jack = FALSE; with the default
# jack = TRUE the interval is a delete-one-donor jackknife and n_boot is not
# consulted at all. Both branches are run, so the supplement can say which
# parameter the interval actually depends on.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(effn))

SM <- Sys.getenv("EFFN_SMALL", "C:/rtmp/effnfig2/gse96583_small.rds")
OUT <- Sys.getenv("EFFN_FIG", "C:/rtmp/effnfig2")
s <- readRDS(SM)
cat("subsample:", nrow(s$counts), "genes x", ncol(s$counts), "cells\n")

NBOOT <- c(10L, 40L, 100L, 200L, 400L)
SEEDS <- c(20261021L, 20261022L, 20261023L)
NTECH <- 40L

rows <- list()
for (jk in c(FALSE, TRUE)) {
  for (b in NBOOT) {
    for (sd in SEEDS) {
      t0 <- proc.time()[["elapsed"]]
      r <- effn(s$counts, group = s$group, donor = s$donor,
                n_boot = b, n_tech = NTECH, seed = sd, jack = jk)
      el <- proc.time()[["elapsed"]] - t0
      sm <- r$summary
      ci <- sm$median_mde_log2_ci
      ok <- r$per_gene[r$per_gene$status == "ok", ]
      rows[[length(rows) + 1L]] <- data.frame(
        jack = jk, n_boot = b, seed = sd,
        median_sigma_tech = stats::median(ok$sigma_tech, na.rm = TRUE),
        median_n_eff = sm$median_n_eff,
        median_mde_fold = 2^sm$median_mde_log2,
        ci_lo = 2^ci[1], ci_hi = 2^ci[2],
        ci_width = ci[2] - ci[1],
        ci_method = sm$ci_method,
        n_ok = sm$n_ok, secs = el)
      cat(sprintf("jack=%-5s n_boot=%3d seed=%d  mde=%.6f  ci=[%.4f, %.4f]  %.1fs\n",
                  jk, b, sd, 2^sm$median_mde_log2, 2^ci[1], 2^ci[2], el))
    }
  }
}

res <- do.call(rbind, rows)
write.table(res, file.path(OUT, "t_nboot_sensitivity.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## per branch: how far the point estimate moves across n_boot, and how much
## the interval ripples. The point estimate should be identical; if it is not,
## the seed streams are coupled and the claim in the manual is wrong.
sm <- do.call(rbind, lapply(split(res, list(res$jack, res$n_boot), drop = TRUE),
  function(d) data.frame(
    jack = d$jack[1], n_boot = d$n_boot[1],
    mde_mean = mean(d$median_mde_fold),
    mde_spread_pct = 100 * (max(d$median_mde_fold) - min(d$median_mde_fold)) /
                       mean(d$median_mde_fold),
    sigma_tech_mean = mean(d$median_sigma_tech),
    n_eff_mean = mean(d$median_n_eff),
    ci_lo_mean = mean(d$ci_lo), ci_hi_mean = mean(d$ci_hi),
    ci_width_mean = mean(d$ci_width),
    ci_width_cv_pct = 100 * sd(d$ci_width) / mean(d$ci_width),
    secs_mean = mean(d$secs))))
sm <- sm[order(sm$jack, sm$n_boot), ]
write.table(sm, file.path(OUT, "t_nboot_sensitivity_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## the number the supplement quotes: the widest interval minus the narrowest,
## as a share of the mean, inside the bootstrap branch
bs <- sm[!sm$jack, ]
cat(sprintf("\nbootstrap branch: ci width %.3f to %.3f log2 (%.1f%% of mean)\n",
            min(bs$ci_width_mean), max(bs$ci_width_mean),
            100 * (max(bs$ci_width_mean) - min(bs$ci_width_mean)) / mean(bs$ci_width_mean)))
cat(sprintf("point estimate across n_boot: max spread %.3g%%\n",
            max(sm$mde_spread_pct)))
jk <- sm[sm$jack, ]
cat(sprintf("jackknife branch: ci width identical across n_boot? %s (%.6f)\n",
            length(unique(round(jk$ci_width_mean, 12))) == 1, jk$ci_width_mean[1]))
cat("written:", file.path(OUT, "t_nboot_sensitivity.tsv"), "\n")
