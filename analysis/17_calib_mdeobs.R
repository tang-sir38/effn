#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 17_calib_mdeobs.R
#
# 诊断一个具体问题：v1.1.2 的校准在 rho 高时 power 冲到 0.95，是哪里来的。
#
# 已知的机制（06_selection_bias.R 已量化）：effn 只对 s2_bio > 0 的基因给 MDE，
# 这等价于要求 s2_obs > s2_tech。rho 越大，被这个条件剔掉的越是"观测方差偏低"
# 的基因，留下来的那批方差系统性偏大，MDE 被抬高。抬高之后，用这个 MDE 去
# 模拟，实际 power 自然超过 0.80。
#
# 但这里有两个效应缠在一起：
#   (a) 选择：进入 MDE 中位数的基因不是随机样本
#   (b) 计算：MDE 用 s2_bio 还是 s2_obs
# 代数上二者对同一个基因恒等（sigma_bio*sqrt(1/n_eff) = sigma_obs/sqrt(n)），
# 所以差别只可能来自 (a)。本脚本把两者分开：
#
#   target 集只按"真值 sigma 分层 + 表达量"选，与估计出来的方差无关 —— 消除 (a)
#   然后对同一批基因分别用 bio 口径和 obs 口径算 MDE，比较
#     ① 各自与真值 MDE 的比值
#     ② 在这个 MDE 上模拟，实际 power 是否回到 0.80
#
# 若 obs 口径把 power 拉回 0.80，说明问题就是选择，修法是"对全部表达基因给 MDE"，
# 而不是"改公式"。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(effn))

OUT <- Sys.getenv("EFFN_OUT", ".")
set.seed(20261021)

G_STRATUM <- 25L
SIGMA_TRUE <- c(0.20, 0.35, 0.50, 0.70)
SIG_HET <- 0.35
G <- G_STRATUM * length(SIGMA_TRUE)
ND <- 8L
TS <- qt(0.975, ND - 1L) + qt(0.80, ND - 1L)

ARMS <- list(
  deep = list(NCELL = 300L, SIZE = 10,  L = 5000, CPM_MAX = 2000, RT = 2000L),
  mid  = list(NCELL = 20L,  SIZE = 2,   L = 800,  CPM_MAX = 500,  RT = 2000L),
  high = list(NCELL = 10L,  SIZE = 1.5, L = 1500, CPM_MAX = 500,  RT = 2000L)
)

run_arm <- function(arm, p) {
  NCELL <- p$NCELL; SIZE_CELL <- p$SIZE; L <- p$L
  CPM <- exp(seq(log(5), log(p$CPM_MAX), length.out = G))
  pvec <- CPM / sum(CPM)
  sig_true <- rep(SIGMA_TRUE, each = G_STRATUM)

  draw_b <- function() {
    list(b = matrix(rnorm(G * ND, 0, sig_true), nrow = G),
         e = matrix(rnorm(G * ND, 0, SIG_HET), nrow = G))
  }
  sim_cells <- function() {
    bb <- draw_b(); idx <- rep(seq_len(ND), each = NCELL)
    yc <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(pvec * exp(bb$b)[, idx] * L), size = SIZE_CELL), nrow = G)
    yt <- matrix(rnbinom(G * ND * NCELL, mu = as.numeric(pvec * exp(bb$b + bb$e)[, idx] * L), size = SIZE_CELL), nrow = G)
    list(counts = cbind(yc, yt), donor = rep(idx, times = 2L),
         group = rep(c("ctrl", "case"), each = ND * NCELL))
  }
  sim_pb <- function(delta_g = rep(0, G)) {
    bb <- draw_b()
    yc <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b) * L * NCELL), size = NCELL * SIZE_CELL), nrow = G)
    yt <- matrix(rnbinom(G * ND, mu = as.numeric(pvec * exp(bb$b + matrix(delta_g, nrow = G, ncol = ND) + bb$e) * L * NCELL), size = NCELL * SIZE_CELL), nrow = G)
    cbind(yc, yt)
  }
  contrast <- function(y) {
    cpm <- log2(sweep(y, 2, colSums(y), "/") * 1e6 + 1)
    cpm[, (ND + 1L):(2L * ND), drop = FALSE] - cpm[, seq_len(ND), drop = FALSE]
  }
  paired_t_p <- function(y) {
    d <- contrast(y); m <- rowMeans(d)
    s <- sqrt(rowSums((d - m)^2) / (ND - 1L))
    tst <- m / (s / sqrt(ND)); 2 * stats::pt(-abs(tst), ND - 1L)
  }

  ## 真值：同一生成过程下对比值的抽样 SD
  RT <- p$RT
  dnull <- array(NA_real_, dim = c(RT * ND, G))
  for (r in seq_len(RT)) dnull[((r - 1L) * ND + 1L):(r * ND), ] <- t(contrast(sim_pb()))
  sd_true <- apply(dnull, 2, stats::sd); rm(dnull)
  mde_true <- TS * sd_true / sqrt(ND)

  one <- sim_cells()
  est <- effn(one$counts, one$group, one$donor, n_boot = 0L, n_tech = 40L, seed = 1L)
  pg <- est$per_gene
  mde_bio <- pg$mde_log2
  mde_obs <- TS * pg$sigma_obs / sqrt(ND)

  ## target 只按真值分层 + 表达量选，不看估计出来的方差
  expressed <- which(!is.na(pg$sigma_obs) & pg$mu_cpm >= 0.1)
  target <- unlist(lapply(SIGMA_TRUE, function(s) {
    i <- expressed[sig_true[expressed] == s]
    i <- i[order(-pg$mu_cpm[i])]
    i[seq_len(min(5L, length(i)))]
  }))

  R <- 300L
  rate <- function(delta) {
    rej <- integer(length(target))
    for (r in seq_len(R)) rej <- rej + (paired_t_p(sim_pb(delta))[target] < 0.05)
    mean(rej / R)
  }
  sim_at <- function(vec) {
    delta <- rep(0, G); delta[target] <- vec[target]
    delta[!is.finite(delta)] <- 0
    delta[-target] <- -sum(pvec[target] * delta[target]) / sum(pvec[-target])
    rate(delta)
  }

  in_ok <- pg$status[target] == "ok"
  cat(sprintf("\n[%s] rho %.3f | targets %d, of which status ok: %d\n",
              arm, est$summary$median_rho_tech, length(target), sum(in_ok)))
  cat(sprintf("[%s] median MDE  true %.3f | bio %.3f | obs %.3f  (log2)\n",
              arm, median(mde_true[target]), median(mde_bio[target], na.rm = TRUE),
              median(mde_obs[target])))
  cat(sprintf("[%s] ratio to truth: bio %.3f | obs %.3f\n", arm,
              median(mde_bio[target] / mde_true[target], na.rm = TRUE),
              median(mde_obs[target] / mde_true[target])))

  pw_bio <- sim_at(mde_bio)
  pw_obs <- sim_at(mde_obs)
  pw_true <- sim_at(mde_true)
  cat(sprintf("[%s] power at own MDE: bio %.3f | obs %.3f | oracle(true) %.3f\n",
              arm, pw_bio, pw_obs, pw_true))

  data.frame(arm = arm, median_rho = est$summary$median_rho_tech,
             n_target = length(target), n_target_ok = sum(in_ok),
             med_mde_true = median(mde_true[target]),
             med_mde_bio = median(mde_bio[target], na.rm = TRUE),
             med_mde_obs = median(mde_obs[target]),
             ratio_bio = median(mde_bio[target] / mde_true[target], na.rm = TRUE),
             ratio_obs = median(mde_obs[target] / mde_true[target]),
             power_bio = pw_bio, power_obs = pw_obs, power_oracle = pw_true)
}

res <- do.call(rbind, lapply(names(ARMS), function(a) run_arm(a, ARMS[[a]])))
print(res, digits = 4)
utils::write.table(res, file.path(OUT, "t_calib_mdeobs.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("\nMDEODONE\n")
