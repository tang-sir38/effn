## 06_selection_bias.R -- 验证四臂校准里 power 偏高的来源
##
## 现象：四臂校准中 power at MDE 随 rho 单调上升
##   deep    rho 0.011 -> power 0.767
##   shallow rho 0.232 -> power 0.874
##   mid     rho 0.527 -> power 0.893
##   high    rho 0.674 -> power 0.948
## 若 MDE 无偏，四条都该在 0.80 附近。
##
## 假设：effn 只对 s2_bio > 0 的基因给 MDE（status == "ok"），而
##       s2_bio = s2_obs - s2_tech > 0 等价于 s2_obs > s2_tech。
##       这是一个以观测方差为条件的选择：rho 越大（s2_tech 越大），
##       被筛掉的越是 s2_obs 偏小的基因，留下的基因 s2_obs 系统性偏大，
##       于是 MDE = t * sqrt(s2_obs) / sqrt(n) 系统性偏大，power 偏高。
##
## 本脚本用同一份模拟数据检验：按 status 分组比较 sigma_obs，
## 并与大样本真值 sd_true 对照。若 ok 组的 sigma_obs / sd_true 中位数
## 显著大于 1，且随 rho 增大而增大，假设成立。

LIB <- "C:/rtmp/effnlib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages(library(effn))

set.seed(20261021)
G_STRATUM <- 25L
SIGMA_TRUE <- c(0.20, 0.35, 0.50, 0.70)
SIG_HET <- 0.35
G <- G_STRATUM * length(SIGMA_TRUE)
ND <- 8L

ARMS <- list(
  deep    = list(NCELL = 300L, SIZE = 10,  L = 5000, CPM_MAX = 2000, RT = 2000L),
  shallow = list(NCELL = 40L,  SIZE = 3,   L = 1500, CPM_MAX = 500,  RT = 1000L),
  mid     = list(NCELL = 20L,  SIZE = 2,   L = 800,  CPM_MAX = 500,  RT = 1000L),
  high    = list(NCELL = 10L,  SIZE = 1.5, L = 1500, CPM_MAX = 500,  RT = 1000L)
)

run <- function(arm, p) {
  NCELL <- p$NCELL; SIZE_CELL <- p$SIZE; L <- p$L
  CPM <- exp(seq(log(5), log(p$CPM_MAX), length.out = G))
  pvec <- CPM / sum(CPM)
  sig_true <- rep(SIGMA_TRUE, each = G_STRATUM)

  draw_b <- function() {
    b <- matrix(rnorm(G * ND, 0, sig_true), nrow = G)
    e <- matrix(rnorm(G * ND, 0, SIG_HET), nrow = G)
    list(b = b, e = e)
  }
  sim_cells <- function() {
    bb <- draw_b()
    idx <- rep(seq_len(ND), each = NCELL)
    yc <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(pvec * exp(bb$b)[, idx] * L), size = SIZE_CELL), nrow = G)
    yt <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(pvec * exp(bb$b + bb$e)[, idx] * L), size = SIZE_CELL), nrow = G)
    list(counts = cbind(yc, yt), donor = rep(idx, times = 2L),
         group = rep(c("ctrl", "case"), each = ND * NCELL))
  }
  sim_pb <- function() {
    bb <- draw_b()
    yc <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b) * L * NCELL), size = NCELL * SIZE_CELL), nrow = G)
    yt <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b + bb$e) * L * NCELL), size = NCELL * SIZE_CELL), nrow = G)
    cbind(yc, yt)
  }
  contrast <- function(y) {
    cpm <- log2(sweep(y, 2, colSums(y), "/") * 1e6 + 1)
    cpm[, (ND + 1L):(2L * ND), drop = FALSE] - cpm[, seq_len(ND), drop = FALSE]
  }

  ## 大样本真值：跨重复的差值 SD
  RT <- p$RT
  dnull <- array(NA_real_, dim = c(RT * ND, G))
  for (r in seq_len(RT)) dnull[((r - 1L) * ND + 1L):(r * ND), ] <- t(contrast(sim_pb()))
  sd_true <- apply(dnull, 2, stats::sd)
  rm(dnull)

  one <- sim_cells()
  est <- effn(one$counts, one$group, one$donor, n_boot = 0L, n_tech = 40L, seed = 1L)
  pg <- est$per_gene
  ok <- pg$status == "ok"

  ## 关键量：sigma_obs 相对于真值的比值，按 status 分组
  r_all <- pg$sigma_obs / sd_true
  data.frame(
    arm = arm, NCELL = NCELL, SIZE = SIZE_CELL, L = L,
    rho = est$summary$median_rho_tech,
    n_ok = sum(ok),
    n_techdom = sum(pg$status == "tech_dominated"),
    n_lowres = sum(pg$status == "low_resolution"),
    n_lowexpr = sum(pg$status == "low_expression"),
    ## 全部基因（无选择）
    ratio_all_median = median(r_all, na.rm = TRUE),
    ## 仅 ok 基因（有条件选择）
    ratio_ok_median = median(r_all[ok], na.rm = TRUE),
    ## 被 tech_dominated 筛掉的基因（应当偏小）
    ratio_techdom_median = median(r_all[pg$status == "tech_dominated"], na.rm = TRUE),
    ## 选择的强度：ok 组相对全部基因的方差膨胀
    var_inflation = median(pg$sigma_obs[ok]^2, na.rm = TRUE) /
                    median(pg$sigma_obs^2, na.rm = TRUE)
  )
}

res <- do.call(rbind, lapply(names(ARMS), function(a) run(a, ARMS[[a]])))
print(res, digits = 4, row.names = FALSE)
write.table(res, "t_selection_bias.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n判读：\n")
cat("ratio_all_median 若 ≈ 1，说明方差估计本身无偏，问题只在选择。\n")
cat("ratio_ok_median 若 > ratio_all_median，说明 ok 筛选抬高了观测方差。\n")
cat("var_inflation > 1 的量即为 MDE 被抬高的倍数（开方后）。\n")
cat("SELBIAS_DONE\n")
