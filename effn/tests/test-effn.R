## Base-R test suite (no testthat dependency).
## Every block ends with a printed "ok" line; any failed stopifnot aborts the file
## and R CMD check reports an ERROR.

suppressPackageStartupMessages(library(effn))

sim_counts <- function(G = 120, n_per = 4, cells = 25, delta = 0.5, seed = 1) {
  set.seed(seed)
  donors <- paste0("d", seq_len(2 * n_per))
  grp <- rep(c("ctrl", "case"), each = n_per)
  cell_donor <- rep(donors, each = cells)
  cell_grp <- rep(grp, each = cells)
  mu <- exp(rnorm(G, mean = 2.2, sd = 0.8))
  b <- matrix(rnorm(G * length(donors), 0, 0.45), nrow = G)
  eff <- ifelse(grp == "case", delta, 0)
  eta <- outer(log(mu), rep(1, length(donors))) + b +
    outer(rep(1, G), eff)
  base <- exp(eta)
  cnt <- matrix(0, nrow = G, ncol = length(cell_donor))
  for (j in seq_along(cell_donor)) {
    cnt[, j] <- rnbinom(G, mu = base[, match(cell_donor[j], donors)], size = 10)
  }
  rownames(cnt) <- paste0("g", seq_len(G))
  colnames(cnt) <- paste0("c", seq_along(cell_donor))
  list(counts = cnt, donor = cell_donor, group = cell_grp)
}

## 1. normal case -----------------------------------------------------------
d <- sim_counts()
r <- effn(d$counts, d$group, d$donor, n_boot = 25, n_tech = 20)
stopifnot(inherits(r, "effn_report"))
stopifnot(nrow(r$per_gene) == nrow(d$counts))
stopifnot(all(c("sigma_bio", "rho_tech", "mde_log2", "status") %in% names(r$per_gene)))
stopifnot(r$summary$n_donor == 8L)
stopifnot(r$summary$n_ok > 20L)
ok <- r$per_gene$status == "ok"
stopifnot(all(r$per_gene$n_eff_group1[ok] <= 4 + 1e-8))
stopifnot(all(r$per_gene$mde_log2[ok] > 0))
## the resolution floor keeps every reported MDE finite, and no estimable gene
## sits below two effective donors
stopifnot(all(is.finite(r$per_gene$mde_log2[ok])))
stopifnot(!any(is.infinite(r$per_gene$mde_log2)))
stopifnot(all(r$per_gene$n_eff[ok] >= 2 - 1e-8))
stopifnot(r$summary$n_unestimable[["low_resolution"]] >= 0L)
cat("ok - normal case\n")

## 2. determinism -----------------------------------------------------------
r2 <- effn(d$counts, d$group, d$donor, n_boot = 25, n_tech = 20)
stopifnot(identical(r$per_gene, r2$per_gene))
stopifnot(identical(r$summary$median_n_eff, r2$summary$median_n_eff))
cat("ok - byte-identical on rerun\n")

## 3. sparse Matrix path agrees with dense path ----------------------------
if (requireNamespace("Matrix", quietly = TRUE)) {
  sm <- Matrix::Matrix(d$counts, sparse = TRUE)
  r3 <- effn(sm, d$group, d$donor, n_boot = 25, n_tech = 20)
  stopifnot(max(abs(r3$per_gene$sigma_bio - r$per_gene$sigma_bio)) < 1e-10)
  cat("ok - sparse and dense agree\n")
}

## 4. all-zero gene is reported, not dropped -------------------------------
cc <- d$counts
cc[1, ] <- 0
r4 <- effn(cc, d$group, d$donor, n_boot = 10, n_tech = 10)
stopifnot(sum(r4$per_gene$status == "all_zero") == 1L)
stopifnot(r4$summary$n_unestimable[["all_zero"]] == 1L)
cat("ok - all-zero gene flagged\n")

## 5. non-integer counts are tolerated -------------------------------------
cn <- d$counts
cn[, 1] <- cn[, 1] + 0.5
r5 <- effn(cn, d$group, d$donor, n_boot = 10, n_tech = 10)
stopifnot(inherits(r5, "effn_report"))
cat("ok - non-integer counts tolerated\n")

## 6. prescribed donors are monotone in the target ------------------------
p1 <- prescribe(r, target_mde = 1)
p2 <- prescribe(r, target_mde = 0.5)
stopifnot(p1$shortfall >= 0)
stopifnot(is.na(p2$n_required) || p2$n_required >= p1$n_required)
stopifnot(is.character(p1$sentence), nchar(p1$sentence) > 20)
cat("ok - prescription monotone\n")

## 7. two donors per group is refused --------------------------------------
d2 <- sim_counts(n_per = 2)
err <- try(effn(d2$counts, d2$group, d2$donor, n_boot = 5, n_tech = 10),
           silent = TRUE)
stopifnot(inherits(err, "try-error"))
cat("ok - n = 2 per group refused\n")

## 8. single group is refused ---------------------------------------------
err <- try(effn(d$counts, rep("a", ncol(d$counts)), d$donor,
                n_boot = 5, n_tech = 10), silent = TRUE)
stopifnot(inherits(err, "try-error"))
cat("ok - single group refused\n")

## 9. missing donor labels are refused, never silently n = 1 ---------------
err <- try(effn(d$counts, d$group, donor = NULL, n_boot = 5, n_tech = 10),
           silent = TRUE)
stopifnot(inherits(err, "try-error"))
stopifnot(grepl("donor labels unavailable", as.character(err)))
cat("ok - missing donor labels refused\n")

## 10. paired design is auto-detected, and refused when forced unpaired ----
## build a paired dataset: every donor contributes ctrl and case cells
set.seed(7)
G <- 80; nd <- 6; cells <- 20
donors <- paste0("p", seq_len(nd))
cd <- rep(rep(donors, each = cells), 2)
cg <- rep(c("ctrl", "case"), each = nd * cells)
mu <- exp(rnorm(G, 2, 0.6))
b <- matrix(rnorm(G * nd, 0, 0.5), G)
base <- exp(outer(log(mu), rep(1, nd)) + b)
cnt <- matrix(0, G, length(cd))
for (j in seq_along(cd)) {
  k <- if (cg[j] == "case") 2 else 1
  cnt[, j] <- rnbinom(G, mu = base[, match(cd[j], donors)] * k, size = 10)
}
rownames(cnt) <- paste0("g", seq_len(G))
rp <- effn(cnt, cg, cd, n_boot = 15, n_tech = 20)
stopifnot(identical(rp$summary$design, "paired"))
stopifnot(rp$summary$n_donor == nd)
stopifnot(all(rp$per_gene$n_eff[rp$per_gene$status == "ok"] <= nd + 1e-8))
err <- try(effn(cnt, cg, cd, n_boot = 5, n_tech = 10, paired = FALSE), silent = TRUE)
stopifnot(inherits(err, "try-error"))
pr <- prescribe(rp, target_mde = log2(1.5))
stopifnot(isTRUE(pr$paired))
cat("ok - paired design detected and handled\n")

## 11. plot returns panel data -------------------------------------------
pdf(NULL)
pd <- plot_effn(r, target_mde = 1)
dev.off()
stopifnot(all(c("neff", "mde", "loss") %in% names(pd)))
cat("ok - plot_effn returns panel data\n")

## 12. the donor solver matches a brute-force scan (no oscillation) -------
for (pr in c(TRUE, FALSE)) {
  for (tg in c(1, log2(1.5), 0.25)) {
    got <- as.integer(effn:::prescribe_n(tg, 0.34, 1.47, paired = pr))
    bf <- NA_integer_
    for (nn in 3:5000) {
      ne <- nn / (1 + 1.47)
      ## design df, not effective df: n - 1 paired, 2n - 2 unpaired
      df <- if (pr) nn - 1 else 2 * nn - 2
      m <- if (df <= 0) Inf else
        (stats::qt(0.975, df) + stats::qt(0.8, df)) * 0.34 *
        (if (pr) sqrt(1 / ne) else sqrt(2 / ne))
      if (m <= tg) { bf <- as.integer(nn); break }
    }
    stopifnot(identical(got, bf))
  }
}
stopifnot(is.na(effn:::prescribe_n(0.01, 0.34, 1.47, paired = TRUE)))
cat("ok - donor solver matches brute force\n")

## 13. the summary and prescribe() agree, and neither says unreachable -----
stopifnot(!is.na(r$summary$n_required_2fold))
stopifnot(identical(as.numeric(p1$n_required), as.numeric(r$summary$n_required_2fold)))
stopifnot(p1$n_required > 0)
stopifnot(!grepl("not reachable", p1$sentence))
cat("ok - summary prescription consistent\n")

## 14. SingleCellExperiment route, when available -------------------------
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = d$counts),
    colData = data.frame(donor = d$donor, group = d$group))
  ## the same arguments must give the same answer through both routes; the
  ## point estimate depends on n_tech, so the comparison has to fix it
  r5 <- effn(d$counts, d$group, d$donor, n_boot = 10, n_tech = 10)
  r6 <- effn(sce, d$group, n_boot = 10, n_tech = 10)
  stopifnot(identical(r5$summary$median_sigma_bio, r6$summary$median_sigma_bio))
  stopifnot(identical(r5$per_gene$n_eff, r6$per_gene$n_eff))
  ## donor read from colData; group given as the column name rather than a vector
  r7 <- effn(sce, "group", n_boot = 10, n_tech = 10)
  stopifnot(identical(r7$summary$median_sigma_bio, r5$summary$median_sigma_bio))
  stopifnot(grepl("^metadata:", r7$summary$donor_source))
  cat("ok - SingleCellExperiment route, vector and column name agree\n")
}

## 15. the 800-cell subsample that ships in inst/extdata ------------------
## Point estimates depend on n_tech only; n_boot moves the confidence
## intervals. The medians below are the ones quoted in the paper.
rds <- system.file("extdata", "gse96583_small.rds", package = "effn")
if (!nzchar(rds)) {
  cat("ok - extdata absent, subsample test skipped\n")
} else {
  small <- readRDS(rds)
  stopifnot(identical(dim(small$counts), c(35635L, 800L)))
  tb <- table(small$donor, small$group)
  stopifnot(ncol(tb) == 2L, all(tb == 50L))
  rs <- effn(small$counts, small$group, small$donor,
             n_boot = 30L, n_tech = 40L, seed = 20261021L)
  stopifnot(abs(rs$summary$median_n_eff - 5.2788858444) < 1e-3)
  stopifnot(abs(rs$summary$median_mde_log2 - 2.7603529065) < 1e-3)
  stopifnot(abs(rs$summary$median_rho_tech - 0.5154713955) < 1e-3)
  ## the reported interval must cover its own point estimate
  stopifnot(rs$summary$median_mde_log2_ci[1] < rs$summary$median_mde_log2,
            rs$summary$median_mde_log2_ci[2] > rs$summary$median_mde_log2)
  stopifnot(rs$summary$median_n_eff_ci[1] < rs$summary$median_n_eff,
            rs$summary$median_n_eff_ci[2] > rs$summary$median_n_eff)
  stopifnot(rs$summary$median_n_eff <= length(unique(small$donor)))
  cat("ok - shipped subsample reproduces its medians\n")
}

## 16. the t quantile uses the design df, not n_eff - 1 -------------------
## Regression test for a defect found in review. mde_t() used df = n_eff - 1,
## which shrinks the degrees of freedom by the factor 1 + rho and inflates the
## MDE in exactly the regime the package exists to describe. With 8 paired
## donors the quantile must be evaluated at 7, whatever rho turns out to be.
set.seed(11)
G <- 120
d16 <- rep(paste0("d", 1:8), each = 40)
g16 <- rep(rep(c("ctrl", "case"), each = 20), 8)   # every donor in both groups
mu16 <- exp(rnorm(G, 3, 0.5))
b16 <- matrix(rnorm(G * 8, 0, 0.4), G)
cnt16 <- sapply(seq_along(d16), function(j)
  rnbinom(G, mu = mu16 * exp(b16[, match(d16[j], paste0("d", 1:8))]), size = 20))
rownames(cnt16) <- paste0("g", 1:G)
r16 <- effn(cnt16, g16, d16, n_boot = 0L, n_tech = 20L, seed = 1)
pg16 <- r16$per_gene[r16$per_gene$status == "ok", ]
manual <- (stats::qt(0.975, 7) + stats::qt(0.8, 7)) * pg16$sigma_bio *
  sqrt(1 / pg16$n_eff)
stopifnot(max(abs(manual - pg16$mde_log2)) < 1e-8)
## the standard error of the contrast is sigma_bio * sqrt(1/n_eff) = sigma_obs /
## sqrt(n), so rho enters only through the quantile
stopifnot(max(abs(pg16$sigma_bio * sqrt(1 / pg16$n_eff) -
                    pg16$sigma_obs / sqrt(8))) < 1e-8)
cat("ok - t quantile uses the design df (n - 1)\n")

cat("ALL TESTS PASSED\n")
