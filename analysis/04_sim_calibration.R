#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 04_sim_calibration.R
# End-to-end calibration. One cell-level dataset is used to estimate; many fresh
# datasets from the same truth test whether the predicted power is realised.
#
# Design: paired, 8 donors, NCELL cells per donor-group, genes spanning a CPM
# grid. Per cell the count is
#   NB(mean = CPM/1e6 * L * exp(b_dg + delta_g * case + eps_dg * case), size = S)
# with b_dg ~ N(0, sigma_bio) the donor level and eps_dg ~ N(0, SIG_HET) the
# donor-specific response. The response term matters: with a single shared
# effect the within-donor contrast has no biological variance left, sigma_bio^2
# collapses onto the technical term, and every gene is reported as
# tech-dominated. That is a property of the generator, not of the estimator.
#
# Two simulators share the same truth: sim_cells() draws every cell and is used
# once, because effn() needs cells to measure the technical variance; sim_pb()
# draws only the donor-group sums as NB(mean = NCELL * ..., size = NCELL * S),
# which is the exact distribution of a sum of NCELL cells of NB(size = S), and
# is used for the replicates because the paired t-test needs nothing finer.
#
# Two arms, because a calibration is only evidence for the regime it visits.
#   deep    300 cells per donor-group, size 10  -> rho near 0.01
#   shallow  40 cells per donor-group, size 3   -> rho of order one, the regime
#           of the worked dataset (rho = 1.05 there)
# The deep arm alone cannot test the depth term: with almost no technical
# variance the depth penalty is a rounding error, so an estimator that inflated
# it by the factor 1 + rho would still pass.
#
# What is checked, per arm:
#   1. sigma_bio matches the empirical SD of the within-donor contrast;
#   2. the realised rejection rate at delta = MDE is 0.80 +/- 0.10;
#   3. the realised rejection rate at delta = 0 is 0.05 +/- 0.02.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(effn))

OUT <- Sys.getenv("EFFN_OUT", ".")
set.seed(20261021)

G_STRATUM <- 25L
SIGMA_TRUE <- c(0.20, 0.35, 0.50, 0.70)
SIG_HET <- 0.35
G <- G_STRATUM * length(SIGMA_TRUE)
ND <- 8L

ARMS <- list(
  deep = list(NCELL = 300L, SIZE = 10, L = 5000, CPM_MAX = 2000, RT = 2000L),
  shallow = list(NCELL = 40L, SIZE = 3, L = 1500, CPM_MAX = 500, RT = 1000L)
)

run_arm <- function(arm, p) {
  NCELL <- p$NCELL; SIZE_CELL <- p$SIZE; L <- p$L
  CPM <- exp(seq(log(5), log(p$CPM_MAX), length.out = G))
  ## the grid has to carry the whole library: without this the genes hold 1% of
  ## each cell's counts and every gene becomes sparse regardless of its CPM
  pvec <- CPM / sum(CPM)
  sig_true <- rep(SIGMA_TRUE, each = G_STRATUM)

  draw_b <- function() {
    b <- matrix(rnorm(G * ND, 0, sig_true), nrow = G)
    e <- matrix(rnorm(G * ND, 0, SIG_HET), nrow = G)
    list(b = b, e = e)
  }

  sim_cells <- function(delta_g = rep(0, G)) {
    bb <- draw_b()
    mu_c <- pvec * exp(bb$b)
    mu_t <- pvec * exp(bb$b + matrix(delta_g, nrow = G, ncol = ND) + bb$e)
    idx <- rep(seq_len(ND), each = NCELL)
    yc <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(mu_c[, idx] * L), size = SIZE_CELL),
                 nrow = G)
    yt <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(mu_t[, idx] * L), size = SIZE_CELL),
                 nrow = G)
    list(counts = cbind(yc, yt),
         donor = rep(idx, times = 2L),
         group = rep(c("ctrl", "case"), each = ND * NCELL))
  }

  sim_pb <- function(delta_g = rep(0, G)) {
    bb <- draw_b()
    yc <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b) * L * NCELL),
                         size = NCELL * SIZE_CELL), nrow = G)
    yt <- matrix(rnbinom(G * ND,
                         mu = as.numeric(pvec * exp(bb$b +
                                                    matrix(delta_g, nrow = G, ncol = ND) +
                                                    bb$e) * L * NCELL),
                         size = NCELL * SIZE_CELL), nrow = G)
    cbind(yc, yt)
  }

  contrast <- function(y) {
    cpm <- log2(sweep(y, 2, colSums(y), "/") * 1e6 + 1)
    cpm[, (ND + 1L):(2L * ND), drop = FALSE] - cpm[, seq_len(ND), drop = FALSE]
  }

  paired_t <- function(y) {
    d <- contrast(y)
    m <- rowMeans(d)
    s <- sqrt(rowSums((d - m)^2) / (ND - 1L))
    tst <- m / (s / sqrt(ND))
    list(p = 2 * stats::pt(-abs(tst), ND - 1L), diff = m)
  }

  ## ---- empirical truth: SD of the within-donor contrast at delta = 0 ------
  RT <- p$RT
  dnull <- array(NA_real_, dim = c(RT * ND, G))
  for (r in seq_len(RT)) {
    dnull[((r - 1L) * ND + 1L):(r * ND), ] <- t(contrast(sim_pb()))
  }
  sd_true <- apply(dnull, 2, stats::sd)
  rm(dnull)
  cat(sprintf("[%s] empirical SD of the contrast: %.3f to %.3f\n",
              arm, min(sd_true), max(sd_true)))

  ## ---- estimate once -----------------------------------------------------
  one <- sim_cells()
  est <- effn(one$counts, one$group, one$donor, n_boot = 0L, n_tech = 40L, seed = 1L)
  cat(sprintf("[%s] estimator: n_ok %d of %d genes\n", arm, est$summary$n_ok, G))
  pg <- est$per_gene
  ok <- pg$status == "ok"
  if (sum(ok) < 10L) stop("too few estimable genes to calibrate")
  ratio <- pg$sigma_bio[ok] / sd_true[ok]
  cat(sprintf("[%s] sigma_bio / empirical SD: median %.3f (range %.3f-%.3f)\n",
              arm, median(ratio), min(ratio), max(ratio)))
  cat(sprintf("[%s] median n_eff %.2f of %d; median rho %.3f; median MDE %.3f log2\n",
              arm, est$summary$median_n_eff, ND, est$summary$median_rho_tech,
              est$summary$median_mde_log2))

  ## ---- realised power ----------------------------------------------------
  R <- 200L
  mde <- pg$mde_log2
  target <- unlist(lapply(SIGMA_TRUE, function(s) {
    i <- which(ok & sig_true == s)
    i[seq_len(min(5L, length(i)))]
  }))
  res <- list()
  for (f in c(0, 0.5, 1.0, 1.5)) {
    delta_g <- rep(0, G)
    delta_g[target] <- f * mde[target]
    delta_g[!is.finite(delta_g)] <- 0
    ## keep the library size fixed: sum_g p_g * delta_g = 0
    delta_g[-target] <- -sum(pvec[target] * delta_g[target]) / sum(pvec[-target])
    rej <- integer(length(target))
    for (r in seq_len(R)) {
      tt <- paired_t(sim_pb(delta_g))
      rej <- rej + (tt$p[target] < 0.05)
    }
    rate <- rej / R
    by_s <- vapply(SIGMA_TRUE, function(s) {
      k <- sig_true[target] == s
      if (sum(k) == 0L) NA_real_ else mean(rate[k])
    }, numeric(1))
    res[[length(res) + 1L]] <- data.frame(
      arm = arm, factor = f, genes = length(target), mean_rejection = mean(rate),
      q05 = quantile(rate, 0.05), q95 = quantile(rate, 0.95),
      by_stratum = paste(sprintf("%.2f:%s", SIGMA_TRUE,
                                 ifelse(is.na(by_s), "na", sprintf("%.3f", by_s))),
                         collapse = " ")
    )
    cat(sprintf("[%s] delta = %.1f x MDE -> mean rejection %.3f\n",
                arm, f, mean(rate)))
  }
  out <- do.call(rbind, res)
  attr(out, "arm_summary") <- data.frame(
    arm = arm, cells_per_donor = NCELL, size = SIZE_CELL,
    median_rho = est$summary$median_rho_tech,
    median_n_eff = est$summary$median_n_eff,
    median_mde_log2 = est$summary$median_mde_log2,
    sigma_ratio_median = median(ratio),
    ## what the decomposition predicts for sigma_bio / sd_true once a fraction
    ## rho/(1+rho) of the contrast variance is technical: 1 / sqrt(1 + rho)
    sigma_ratio_expected = 1 / sqrt(1 + est$summary$median_rho_tech),
    power_at_mde = out$mean_rejection[out$factor == 1],
    type1_at_zero = out$mean_rejection[out$factor == 0],
    n_ok = sum(ok)
  )
  list(out = out, est = est, ratio = ratio, sd_true = sd_true, pg = pg)
}

runs <- lapply(names(ARMS), function(a) run_arm(a, ARMS[[a]]))
names(runs) <- names(ARMS)

tab <- do.call(rbind, lapply(runs, `[[`, "out"))
arms <- do.call(rbind, lapply(runs, function(r) attr(r$out, "arm_summary")))
rownames(arms) <- NULL

print(tab[, c("arm", "factor", "genes", "mean_rejection", "q05", "q95")])
print(arms, digits = 4)

ok_power <- all(abs(arms$power_at_mde - 0.80) <= 0.10)
ok_type1 <- all(abs(arms$type1_at_zero - 0.05) <= 0.02)
## sigma_bio is compared with the value the decomposition implies, not with the
## raw contrast SD: once 1 - 1/sqrt(1+rho) of that SD is technical, a ratio of
## 1 would mean the estimator had absorbed the technical term into biology.
ok_sig <- all(abs(arms$sigma_ratio_median - arms$sigma_ratio_expected) <= 0.10)
cat(sprintf("\npower at MDE   : %s (target 0.80, tolerance 0.10)\n",
            if (ok_power) "PASS" else "FAIL"))
cat(sprintf("type I at zero : %s (target 0.05, tolerance 0.02)\n",
            if (ok_type1) "PASS" else "FAIL"))
cat(sprintf("sigma_bio      : %s (against 1/sqrt(1+rho) per arm, tolerance 0.10)\n",
            if (ok_sig) "PASS" else "FAIL"))

utils::write.table(tab, file.path(OUT, "t_sim_calibration.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
utils::write.table(arms, file.path(OUT, "t_calibration_arms.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
saveRDS(lapply(runs, function(r) list(summary = r$est$summary, ratio = r$ratio,
                                      per_gene = r$pg)),
        file.path(OUT, "sim_calibration.rds"))
cat("CALIBDONE\n")
