## 15_limma_sim.R -- limma 对照（有真值版）
##
## 08/14 的 limma 对照要等 GSE96583 原始矩阵下载完。本脚本先用模拟数据做一遍，
## 好处是**知道真值**：可以判断 effn 与 limma 谁更接近真实的 MDE，而不只是
## 比较两者差多少。
##
## 三个估计器，同一份数据：
##   effn      供体级 pseudobulk，普通 t，df = n-1
##   limma-PB  供体级 pseudobulk，moderated t，df = df.total
##   limma-cell 细胞级 voom，moderated t
## 参照真值：
##   MDE_true = (t_{1-a/2, n-1} + t_{power, n-1}) * sd_true / sqrt(n)
## 其中 sd_true 由大样本模拟得到的跨重复差值 SD（不是单次样本的 s）。
##
## 预期：(a) limma-PB 的 MDE 小于 effn（收缩 + 先验自由度），但幅度有限；
##       (b) limma-cell 的 MDE 远小于两者，因为它把细胞当重复；
##       (c) effn 与 MDE_true 最接近——这正是它的卖点。
## 若 (a) 的幅度很大（例如 > 20%），说明论文必须把 limma 对照写进正文。

LIB <- "C:/rtmp/effnlib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(limma); library(edgeR); library(effn)
})

set.seed(20261028)
G_STRATUM <- 25L
SIGMA_TRUE <- c(0.20, 0.35, 0.50, 0.70)
SIG_HET <- 0.35
G <- G_STRATUM * length(SIGMA_TRUE)
ND <- 8L
ALPHA <- 0.05; POWER <- 0.80

P <- list(NCELL = 20L, SIZE = 2, L = 800, CPM_MAX = 500)   # mid 臂，rho ~ 0.5

CPM <- exp(seq(log(5), log(P$CPM_MAX), length.out = G))
pvec <- CPM / sum(CPM)
sig_true <- rep(SIGMA_TRUE, each = G_STRATUM)

draw_b <- function() {
  list(b = matrix(rnorm(G * ND, 0, sig_true), nrow = G),
       e = matrix(rnorm(G * ND, 0, SIG_HET), nrow = G))
}
sim_cells <- function() {
  bb <- draw_b(); idx <- rep(seq_len(ND), each = P$NCELL)
  yc <- matrix(rnbinom(G * ND * P$NCELL, mu = as.numeric(pvec * exp(bb$b)[, idx] * P$L), size = P$SIZE), nrow = G)
  yt <- matrix(rnbinom(G * ND * P$NCELL, mu = as.numeric(pvec * exp(bb$b + bb$e)[, idx] * P$L), size = P$SIZE), nrow = G)
  m <- cbind(yc, yt)
  ## 基因名必须在源头就加上：effn 的 per_gene$gene 取自 rownames(counts)，
  ## 若此处为 NULL，后面三路结果的基因名无法对齐。
  rownames(m) <- sprintf("g%04d", seq_len(G))
  list(counts = m, donor = rep(idx, times = 2L),
       group = rep(c("ctrl", "case"), each = ND * P$NCELL))
}
sim_pb <- function() {
  bb <- draw_b()
  yc <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b) * P$L * P$NCELL), size = P$NCELL * P$SIZE), nrow = G)
  yt <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b + bb$e) * P$L * P$NCELL), size = P$NCELL * P$SIZE), nrow = G)
  cbind(yc, yt)
}
contrast <- function(y) {
  cpm <- log2(sweep(y, 2, colSums(y), "/") * 1e6 + 1)
  cpm[, (ND + 1L):(2L * ND), drop = FALSE] - cpm[, seq_len(ND), drop = FALSE]
}

## ---- 真值：跨 1000 次重复的差值 SD ----
RT <- 1000L
dnull <- array(NA_real_, dim = c(RT * ND, G))
for (r in seq_len(RT)) dnull[((r - 1L) * ND + 1L):(r * ND), ] <- t(contrast(sim_pb()))
sd_true <- apply(dnull, 2, stats::sd); rm(dnull)
t_sum <- function(df) qt(1 - ALPHA / 2, df) + qt(POWER, df)
mde_true <- t_sum(ND - 1L) * sd_true / sqrt(ND)
names(mde_true) <- sprintf("g%04d", seq_len(G))   # 与 cnt 的 rownames 对齐
cat(sprintf("true MDE (log2): median %.3f  (%.2f fold)\n", median(mde_true), 2^median(mde_true)))

## ---- 观测数据：一次抽样，三种估计 ----
d <- sim_cells()
cnt <- as.matrix(d$counts)
## 必须给基因名：effn 的 per_gene$gene、limma 的 rownames(pb)、mde_cell 的名字
## 三处要能对上，否则 Reduce(intersect) 得到空集（第一次跑就是这样失败的）。
rownames(cnt) <- sprintf("g%04d", seq_len(nrow(cnt)))
donor <- factor(d$donor); grp <- factor(d$group, levels = c("ctrl", "case"))

## (a) effn
est <- effn(d$counts, d$group, d$donor, n_boot = 0L, n_tech = 40L, seed = 1L)
pg <- est$per_gene
ok <- pg$status == "ok"
cat("effn ok genes:", sum(ok), " median rho:", round(est$summary$median_rho_tech, 3), "\n")

## (b) limma on pseudobulk
key <- paste(donor, grp, sep = "\r"); ug <- unique(key)
pb <- sapply(ug, function(k) Matrix::rowSums(cnt[, key == k, drop = FALSE]))
kd <- factor(sub("\r.*$", "", ug)); kg <- factor(sub("^.*\r", "", ug), levels = levels(grp))
o <- order(kd, kg); pb <- pb[, o, drop = FALSE]; kd <- droplevels(kd[o]); kg <- droplevels(kg[o])
dge <- DGEList(counts = pb); dge <- calcNormFactors(dge, method = "TMM")
des <- model.matrix(~ kd + kg); v <- voom(dge, des, plot = FALSE)
fit <- eBayes(lmFit(v, des)); co <- ncol(des)
se_pb <- fit$stdev.unscaled[, co] * sqrt(fit$s2.post)
mde_pb <- t_sum(fit$df.total) * se_pb
cat(sprintf("limma-PB: df.prior %.2f  df.total %.1f\n", fit$df.prior, median(fit$df.total)))

## (c) limma at cell level
keep <- filterByExpr(cnt, group = grp, min.count = 10)
dge2 <- DGEList(counts = cnt[keep, , drop = FALSE]); dge2 <- calcNormFactors(dge2, method = "TMM")
des2 <- model.matrix(~ donor + grp); v2 <- voom(dge2, des2, plot = FALSE)
fit2 <- eBayes(lmFit(v2, des2))
se_cell <- fit2$stdev.unscaled[, co] * sqrt(fit2$s2.post)
mde_cell <- t_sum(fit2$df.total) * se_cell
cat(sprintf("limma-cell: kept %d genes, df.prior %.2f, df.total %.1f\n",
            sum(keep), fit2$df.prior, median(fit2$df.total)))

## ---- 合并 ----
g_use <- pg$gene[ok]
g_use <- Reduce(intersect, list(g_use, rownames(pb), names(mde_cell)))
cat("genes in all three:", length(g_use), "\n")

res <- data.frame(
  gene = g_use,
  mde_true = mde_true[match(g_use, pg$gene)],
  mde_effn = pg$mde_log2[match(g_use, pg$gene)],
  sigma_obs_effn = pg$sigma_obs[match(g_use, pg$gene)],
  mde_limma_pb = mde_pb[match(g_use, rownames(pb))],
  mde_limma_cell = mde_cell[match(g_use, names(mde_cell))]
)
res$effn_vs_true  <- res$mde_effn / res$mde_true
res$pb_vs_true    <- res$mde_limma_pb / res$mde_true
res$cell_vs_true  <- res$mde_limma_cell / res$mde_true

cat("\n=== MDE 中位数 ===\n")
tab <- data.frame(
  estimator = c("truth (large-sample)", "effn", "limma pseudobulk", "limma cell-level"),
  median_log2 = c(median(res$mde_true), median(res$mde_effn),
                  median(res$mde_limma_pb), median(res$mde_limma_cell))
)
tab$median_fold <- 2^tab$median_log2
tab$vs_truth <- tab$median_fold / tab$median_fold[1]
print(tab, digits = 4, row.names = FALSE)

cat(sprintf("\neffn / truth      : median %.3f  IQR %.3f-%.3f\n",
            median(res$effn_vs_true), quantile(res$effn_vs_true, .25), quantile(res$effn_vs_true, .75)))
cat(sprintf("limma-PB / truth  : median %.3f  IQR %.3f-%.3f\n",
            median(res$pb_vs_true), quantile(res$pb_vs_true, .25), quantile(res$pb_vs_true, .75)))
cat(sprintf("limma-cell / truth: median %.3f  IQR %.3f-%.3f\n",
            median(res$cell_vs_true), quantile(res$cell_vs_true, .25), quantile(res$cell_vs_true, .75)))

write.table(res, "t_limma_sim.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
sink("limma_sim.log")
cat("effn v1.2.0 -- limma comparison on simulated data with known truth\n")
cat("R:", R.version.string, " limma:", as.character(packageVersion("limma")), "\n")
cat("design: ND =", ND, " NCELL =", P$NCELL, " SIZE =", P$SIZE, " L =", P$L, "\n")
cat("seed 20261028; truth from", RT, "replications\n\n")
print(tab, digits = 6, row.names = FALSE)
cat("\nratios vs truth:\n")
print(summary(res[, c("effn_vs_true", "pb_vs_true", "cell_vs_true")]), digits = 4)
sink()
cat("\nLIMMASIM_DONE\n")
