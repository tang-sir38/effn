# =============================================================================
# 00_definitions.R — 课题 A 数学内核的参考实现（只有定义，不跑数据）
# 对应文档：00_ADMIN/DEFINITIONS.md v1.0
# 本文件是"口径"的可执行版本；包内实现见 05_PACKAGE/_support/code/effn/R/*.R，
# 二者由 tests/testthat/test-parity.R 在固定输入上逐值比对，不允许漂移。
# 依赖：base + stats（不依赖 Matrix；稀疏矩阵请先 as.matrix 或传 dgCMatrix 给
#       effn 包的实现，本参考实现只用于口径核对）
# =============================================================================

## ---- 0. 规范化：先求和再归一化 -------------------------------------------
# counts: genes x cells 数值矩阵；donor: 长度 = ncol(counts)
# 返回 donors x genes 的 log2(CPM+1)
pb_norm <- function(counts, donor) {
  stopifnot(is.matrix(counts) || inherits(counts, "Matrix"))
  stopifnot(length(donor) == ncol(counts))
  donor <- as.character(donor)
  ug <- sort(unique(donor))
  # 按 donor 求和：用行索引构造指示矩阵，避免循环
  ind <- matrix(0, nrow = length(ug), ncol = length(donor),
                dimnames = list(ug, NULL))
  ind[cbind(match(donor, ug), seq_along(donor))] <- 1
  c_dg <- as.matrix(ind %*% t(as.matrix(counts)))          # donors x genes
  L_d  <- rowSums(c_dg)
  if (any(L_d <= 0)) stop("donor with zero library size: ", paste(ug[L_d <= 0], collapse = ","))
  y <- log2(sweep(c_dg, 1, L_d, "/") * 1e6 + 1)
  dimnames(y) <- list(ug, rownames(counts))
  y
}

## ---- 1. 细胞层 CPM 均值（用于 B2 低表达判定）-----------------------------
cell_mean_cpm <- function(counts, chunk = 2000L) {
  counts <- as.matrix(counts)
  L <- colSums(counts)
  if (any(L <= 0)) stop("cell with zero library size")
  # 分块，避免大矩阵一次性复制
  g <- nrow(counts); out <- numeric(g)
  for (s in seq(1L, g, by = chunk)) {
    e <- min(s + chunk - 1L, g)
    out[s:e] <- rowMeans(sweep(counts[s:e, , drop = FALSE], 2, L, "/") * 1e6)
  }
  out
}

## ---- 2. 组内 pooled 方差（df = n1 + n2 - 2）------------------------------
sigma_pooled <- function(y, group) {
  # y: donors x genes；group: 长度 = nrow(y) 的因子/字符，必须恰好 2 组
  group <- as.character(group)
  stopifnot(length(group) == nrow(y), length(unique(group)) == 2)
  gs <- sort(unique(group))
  idx <- lapply(gs, function(k) which(group == k))
  n_k <- vapply(idx, length, integer(1))
  if (any(n_k < 2)) stop("each group needs >= 2 donors to compute a variance")
  v_k <- lapply(idx, function(ii) {
    if (length(ii) < 2) return(rep(NA_real_, ncol(y)))
    apply(y[ii, , drop = FALSE], 2, var)
  })
  df <- sum(n_k - 1)
  s2 <- Reduce(`+`, Map(function(v, n) v * (n - 1), v_k, n_k)) / df
  list(s2 = s2, df = df, n = stats::setNames(n_k, gs))
}

## ---- 3. σ_tech：split-half，无分布假设 -----------------------------------
# 返回 list(t_dg = donors x genes 的 σ_tech²(N_d) 估计, usable = 逻辑向量)
sigma_tech_splithalf <- function(counts, donor, K = 20L, seed = 20261021L) {
  donor <- as.character(donor)
  ug <- sort(unique(donor))
  set.seed(seed)
  t_mat <- matrix(NA_real_, nrow = length(ug), ncol = nrow(counts),
                  dimnames = list(ug, rownames(counts)))
  for (d in ug) {
    ci <- which(donor == d)
    Nd <- length(ci)
    if (Nd < 4L) next                       # 无法对半拆
    acc <- numeric(nrow(counts))
    for (k in seq_len(K)) {
      perm <- sample(ci, Nd)
      h <- floor(Nd / 2)
      A <- perm[seq_len(h)]; B <- perm[(h + 1L):Nd]
      yA <- pb_norm(counts[, A, drop = FALSE], rep("a", length(A)))[1, ]
      yB <- pb_norm(counts[, B, drop = FALSE], rep("b", length(B)))[1, ]
      acc <- acc + (yA - yB)^2 / 4          # 每次拆分 = 1 df 的无偏估计
    }
    t_mat[d, ] <- acc / K
  }
  list(t_dg = t_mat, usable = rowSums(!is.na(t_mat)) > 0)
}

## ---- 4. ρ_tech 与 n_eff --------------------------------------------------
rho_and_neff <- function(s2_obs, s2_tech, n1, n2, eps = 0) {
  s2_bio <- pmax(s2_obs - s2_tech, eps)
  rho <- ifelse(s2_bio > 0, s2_tech / s2_bio, NA_real_)
  neff1 <- n1 / (1 + rho); neff2 <- n2 / (1 + rho)
  if (any(neff1 > n1 + 1e-8, na.rm = TRUE) || any(neff2 > n2 + 1e-8, na.rm = TRUE))
    stop("n_eff exceeds n_donor: rho_tech must be >= 0")
  list(s2_bio = s2_bio, sigma_bio = sqrt(s2_bio), rho = rho,
       neff1 = neff1, neff2 = neff2)
}

## ---- 5. MDE（t 精确；禁用 z）---------------------------------------------
mde_t <- function(sigma_bio, neff1, neff2, alpha = 0.05, power = 0.80) {
  df <- neff1 + neff2 - 2
  tt <- stats::qt(1 - alpha / 2, df) + stats::qt(power, df)
  tt * sigma_bio * sqrt(1 / neff1 + 1 / neff2)
}

## ---- 6. 处方：给定 target_mde 反解每组 donor 数 --------------------------
prescribe_n <- function(target_mde, sigma_bio, rho, n_start = 3L,
                        alpha = 0.05, power = 0.80, nmax = 10000L) {
  if (!is.finite(sigma_bio) || sigma_bio <= 0) return(NA_integer_)
  f <- function(n) {
    ne <- n / (1 + rho)
    mde_t(sigma_bio, ne, ne, alpha, power)
  }
  n <- max(3L, as.integer(n_start))
  for (i in seq_len(200)) {
    m <- f(n)
    if (!is.finite(m) || m <= 0) return(NA_integer_)
    n_new <- ceiling(2 * (1 + rho) * ((stats::qt(1 - alpha/2, 2*n/(1+rho) - 2) +
                                        stats::qt(power, 2*n/(1+rho) - 2)) * sigma_bio / target_mde)^2)
    n_new <- max(3L, min(as.integer(n_new), nmax + 1L))
    if (n_new == n) return(n)
    n <- n_new
    if (n > nmax) return(NA_integer_)
  }
  NA_integer_
}

## ---- 7. 三个报告量 -------------------------------------------------------
loss_donor <- function(n_nom, n_donor) n_nom / n_donor
loss_depth <- function(rho) 1 + rho
loss_depth_share <- function(rho) rho / (1 + rho)

## ---- 8. 自检：定义层的断言（run_all.R 会调用）----------------------------
definitions_selfcheck <- function(verbose = TRUE) {
  # (a) pb_norm 与手算一致
  set.seed(1)
  m <- matrix(rpois(60, 5), nrow = 6,
              dimnames = list(paste0("g", 1:6), paste0("c", 1:10)))
  dn <- rep(c("d1", "d2"), each = 5)
  y <- pb_norm(m, dn)
  c1 <- rowSums(m[, 1:5]); L1 <- sum(c1)
  stopifnot(abs(y["d1", 1] - log2(c1[1] / L1 * 1e6 + 1)) < 1e-12)
  # (b) t 分位 > z（n=3/组）
  stopifnot((qt(0.975, 4) + qt(0.8, 4)) > (qnorm(0.975) + qnorm(0.8)))
  # (c) n_eff 恒 <= n_donor
  r <- rho_and_neff(s2_obs = c(1, 2), s2_tech = c(0.2, 0.5), n1 = 8, n2 = 8)
  stopifnot(all(r$neff1 <= 8), all(r$neff2 <= 8))
  # (d) 处方单调：target 越严 -> n 越大
  n_a <- prescribe_n(1.0, 0.8, 0.2); n_b <- prescribe_n(0.5, 0.8, 0.2)
  stopifnot(is.na(n_b) || n_b >= n_a)
  if (verbose) cat("definitions_selfcheck: OK\n")
  invisible(TRUE)
}
