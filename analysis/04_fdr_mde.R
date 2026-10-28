## 04_fdr_mde.R -- FDR-aware MDE 与所需供体数
##
## 目的：v1.1.2 的主数字（中位 MDE 1.63 倍、2 倍变化需 6 个供体）全部建立在
##       名义 alpha = 0.05 上。真实使用场景要控制 FDR。本脚本把同一套 MDE 算法
##       换成"控制 FDR 的等效 alpha"，量化主数字的变化。
##
## 等效 alpha 的来历：Jung (2005) Bioinformatics 21(20):3976-3982 的样本量公式。
##   设 pi0 = 真无差异基因比例，alpha = 单基因检验水准，power = 1-beta，
##   则 FDR = pi0*alpha / (pi0*alpha + (1-pi0)*power)。
##   要求 FDR <= q 并取等号：
##       alpha_eff = q * (1-pi0) * power / (pi0 * (1-q))
##   这个式子的含义：FDR 是"后验"指标，要做的事前控制必须比 q 更严，
##   严格程度由预期信号比例 pi1 = 1-pi0 决定。
##
## MDE 模型（与 v1.1.2 完全一致，见 tests/test-effn.R 第 16 块）：
##   sigma_d  = sigma_bio * sqrt(1 + rho_tech)      # 差值尺度的标准差
##   MDE(n)   = (t_{1-alpha/2, n-1} + t_{power, n-1}) * sigma_d / sqrt(n)
##   自由度取设计自由度 n-1（配对设计），不取 n_eff-1。
##
## 输入：03_RESULTS/gse96583/gse96583_per_gene.tsv.gz
## 输出：03_RESULTS/gse96583/fdr_mde_summary.tsv
##       03_RESULTS/gse96583/fdr_mde_per_gene.tsv.gz
##       03_RESULTS/gse96583/fdr_mde.log

suppressWarnings({
  d <- read.delim(gzfile("gse96583_per_gene.tsv.gz"), stringsAsFactors = FALSE)
})

ALPHA   <- 0.05
POWER   <- 0.80
Q_FDR   <- 0.05
PI1_GRID <- c(0.10, 0.20, 0.30, 0.50)
N_MIN   <- 3L
N_MAX   <- 500L
N_MAIN  <- 8L          # 本研究实际供体数

ok <- d[d$status == "ok" & is.finite(d$sigma_bio) & is.finite(d$rho_tech), ]
cat("genes with status=ok and finite sigma:", nrow(ok), "\n")

ok$sigma_d <- ok$sigma_bio * sqrt(1 + ok$rho_tech)

## ---- 1. 复现性校验：n = 8 时 MDE 应与表内 mde_log2 一致 ----
coef_at <- function(n, alpha) (qt(1 - alpha / 2, n - 1) + qt(POWER, n - 1)) / sqrt(n)
pred8 <- coef_at(N_MAIN, ALPHA) * ok$sigma_d
dev <- abs(pred8 - ok$mde_log2)
cat(sprintf("reproduce check at n=8: max|dev|=%.3e  median|dev|=%.3e  n>1e-8: %d\n",
            max(dev), median(dev), sum(dev > 1e-8)))

## ---- 2. 所需供体数：给定目标 MDE（log2 尺度），反解最小 n ----
nseq <- N_MIN:N_MAX
nreq_from_coef <- function(coefs, thr) {
  ## coefs 单调递减。要最小 n 使 coefs[n] <= thr（与 v1.1.2 的 bisection 判据一致）。
  ## findInterval(-thr, -coefs, left.open = TRUE) = #{i : -coefs[i] < -thr} = #{i : coefs[i] > thr}
  ## 第一个满足 coefs <= thr 的位置是 k + 1。
  ## 注意：不能写成 nseq[k]（那等价于严格判据 coefs < thr，会每个基因少算 1 个供体）。
  k <- findInterval(-thr, -coefs, left.open = TRUE)
  out <- nseq[pmin(k + 1L, length(nseq))]
  out[thr < min(coefs)] <- NA_integer_    # 即使 n = N_MAX 也不够
  as.integer(out)
}

targets <- c(two_fold = log2(2), one_point_five = log2(1.5))

## ---- 3. pi1 网格上的解析汇总 ----
rows <- list()
for (pi1 in PI1_GRID) {
  pi0 <- 1 - pi1
  a_eff <- Q_FDR * pi1 * POWER / (pi0 * (1 - Q_FDR))
  mult  <- coef_at(N_MAIN, a_eff) / coef_at(N_MAIN, ALPHA)
  mde_fdr <- 2^(median(ok$mde_log2) * mult)

  cf_nom <- vapply(nseq, coef_at, numeric(1), alpha = ALPHA)
  cf_fdr <- vapply(nseq, coef_at, numeric(1), alpha = a_eff)

  r <- data.frame(pi1 = pi1, alpha_eff = a_eff, mde_multiplier = mult,
                  median_mde_log2_nominal = median(ok$mde_log2),
                  median_mde_fold_nominal = 2^median(ok$mde_log2),
                  median_mde_log2_fdr = median(ok$mde_log2) * mult,
                  median_mde_fold_fdr = mde_fdr)

  for (tn in names(targets)) {
    n_nom <- nreq_from_coef(cf_nom, targets[[tn]] / ok$sigma_d)
    n_fdr <- nreq_from_coef(cf_fdr, targets[[tn]] / ok$sigma_d)
    r[[paste0("n_", tn, "_nominal_median")]] <- median(n_nom, na.rm = TRUE)
    r[[paste0("n_", tn, "_nominal_q25")]]    <- quantile(n_nom, .25, na.rm = TRUE, names = FALSE)
    r[[paste0("n_", tn, "_nominal_q75")]]    <- quantile(n_nom, .75, na.rm = TRUE, names = FALSE)
    r[[paste0("n_", tn, "_fdr_median")]]     <- median(n_fdr, na.rm = TRUE)
    r[[paste0("n_", tn, "_fdr_q25")]]        <- quantile(n_fdr, .25, na.rm = TRUE, names = FALSE)
    r[[paste0("n_", tn, "_fdr_q75")]]        <- quantile(n_fdr, .75, na.rm = TRUE, names = FALSE)
    r[[paste0("n_", tn, "_fdr_unreachable")]] <- sum(is.na(n_fdr))
    if (tn == "two_fold") {
      per_gene_nom <- n_nom; per_gene_fdr <- n_fdr
    }
  }
  rows[[length(rows) + 1L]] <- r
}
summ <- do.call(rbind, rows)

## ---- 4. 主设定 pi1 = 0.2 的逐基因表 ----
main <- summ[summ$pi1 == 0.20, ]

pg <- data.frame(gene = ok$gene, mu_cpm = ok$mu_cpm,
                 sigma_bio = ok$sigma_bio, rho_tech = ok$rho_tech,
                 sigma_d = ok$sigma_d,
                 mde_fold_nominal = ok$mde_fold,
                 mde_fold_fdr = 2^(ok$mde_log2 * main$mde_multiplier),
                 n_req_2fold_nominal = per_gene_nom,
                 n_req_2fold_fdr = per_gene_fdr)
pg$mde_exceeds_2fold_fdr <- pg$mde_fold_fdr > 2

write.table(summ, "fdr_mde_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
gz <- gzfile("fdr_mde_per_gene.tsv.gz", "w"); write.table(pg, gz, sep = "\t", quote = FALSE, row.names = FALSE); close(gz)

## ---- 5. 报告 ----
cat("\n=== pi1 grid (q = 0.05, power = 0.80, df = 7) ===\n")
print(summ[, c("pi1", "alpha_eff", "mde_multiplier",
               "median_mde_fold_nominal", "median_mde_fold_fdr",
               "n_two_fold_nominal_median", "n_two_fold_fdr_median",
               "n_one_point_five_nominal_median", "n_one_point_five_fdr_median")],
      row.names = FALSE, digits = 4)

cat("\n=== main scenario pi1 = 0.2 ===\n")
cat(sprintf("alpha_eff            %.5f\n", main$alpha_eff))
cat(sprintf("MDE multiplier       %.4f  (on log2 scale)\n", main$mde_multiplier))
cat(sprintf("median MDE nominal   %.4f fold\n", main$median_mde_fold_nominal))
cat(sprintf("median MDE FDR       %.4f fold\n", main$median_mde_fold_fdr))
cat(sprintf("2-fold donors  nom   %g (IQR %g-%g)\n", main$n_two_fold_nominal_median,
            main$n_two_fold_nominal_q25, main$n_two_fold_nominal_q75))
cat(sprintf("2-fold donors  FDR   %g (IQR %g-%g), unreachable %d\n", main$n_two_fold_fdr_median,
            main$n_two_fold_fdr_q25, main$n_two_fold_fdr_q75, main$n_two_fold_fdr_unreachable))
cat(sprintf("1.5-fold donors nom  %g\n", main$n_one_point_five_nominal_median))
cat(sprintf("1.5-fold donors FDR  %g\n", main$n_one_point_five_fdr_median))
cat(sprintf("genes with FDR-aware MDE > 2 fold: %d / %d (%.1f%%)\n",
            sum(pg$mde_exceeds_2fold_fdr), nrow(pg),
            100 * mean(pg$mde_exceeds_2fold_fdr)))

sink("fdr_mde.log")
cat("effn v1.2.0 -- FDR-aware MDE\n")
cat("R:", R.version.string, "\n")
cat("input: gse96583_per_gene.tsv.gz; n_genes_ok =", nrow(ok), "\n")
cat("alpha =", ALPHA, " power =", POWER, " q =", Q_FDR, "\n")
cat("alpha_eff formula: Jung 2005, q*pi1*power/((1-pi1)*(1-q))\n")
cat("reproduce max|dev| at n=8:", max(dev), "\n\n")
print(summ, row.names = FALSE, digits = 6)
sink()
cat("\nwritten: fdr_mde_summary.tsv, fdr_mde_per_gene.tsv.gz, fdr_mde.log\n")
