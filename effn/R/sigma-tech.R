# Technical variance of the donor pseudobulk mean.
#
# Two estimators are provided. `bootstrap` (the default) resamples the donor's
# own cells with replacement at the observed depth and takes the variance of the
# resulting pseudobulk values. It measures Var(pseudobulk mean) directly, on
# whatever scale the pseudobulk uses, and therefore needs no scaling assumption.
#
# `splithalf` divides the cells into two random halves and uses (yA - yB)^2 / 4.
# That expression estimates the variance of a pseudobulk mean at half the
# observed depth; carrying it to the observed depth assumes Var falls as 1/n.
# The assumption holds on the raw count or CPM scale but not on log2(CPM + 1),
# where the transform saturates for sparse genes. Measured on GSE96583, the
# technical SD on log2(CPM + 1) falls 3.0-fold between 25 and 1,000 cells per
# donor (2.03 to 0.68), where a 1/n variance law predicts 6.3-fold, so the
# 1/n step imports a depth trend into rho. The estimator is kept because it is
# cheap and because agreement between the two on the CPM scale is a useful check.
#
# Implementation note: all B resamples of one donor-group are encoded as a
# single sparse B-column weight matrix and applied in one product on that
# group's columns. A dense indicator would force the sparse matrix onto a dense
# path (35,635 x 24,679 would need ~7 GB per call), and a full-length indicator
# per resample re-reads every other group's columns.

#' @keywords internal
#' @noRd
sigma_tech_boot <- function(counts, donor, B = 40L, seed = 20261021L) {
  donor <- as.character(donor)
  ug <- unique(donor)
  set.seed(seed)
  G <- nrow(counts)
  B <- as.integer(B)
  t_mat <- matrix(NA_real_, nrow = length(ug), ncol = G,
                  dimnames = list(ug, rownames(counts)))
  by_group <- split(seq_len(ncol(counts)), donor)

  for (d in ug) {
    ci <- by_group[[d]]
    Nd <- length(ci)
    if (Nd < 4L) next
    ## B multinomial resamples of Nd cells each: counts per cell, columns sum to Nd
    w <- stats::rmultinom(B, size = Nd, prob = rep(1 / Nd, Nd))
    W <- Matrix::sparseMatrix(i = rep(seq_len(Nd), times = B),
                              j = rep(seq_len(B), each = Nd),
                              x = as.numeric(w),
                              dims = c(Nd, B))
    s <- as.matrix(counts[, ci, drop = FALSE] %*% W)
    tot <- colSums(s)
    if (any(tot <= 0)) next
    y <- log2(sweep(s, 2, tot, "/") * 1e6 + 1)
    ## variance across the B pseudobulk values, one per gene
    mu <- rowMeans(y)
    t_mat[d, ] <- rowSums((y - mu)^2) / (B - 1L)
  }
  list(t_dg = t_mat, usable = rowSums(!is.na(t_mat)) > 0)
}

#' Technical variance of the donor pseudobulk mean
#'
#' @param counts genes x cells matrix or Matrix.
#' @param donor donor vector, same length as ncol(counts).
#' @param K number of resamples (bootstrap) or random splits (split-half).
#' @param seed integer seed; the same seed always gives the same answer.
#' @return list with `t_dg` (donors x genes) and `usable` (logical).
#' @keywords internal
#' @noRd
sigma_tech_splithalf <- function(counts, donor, K = 20L, seed = 20261021L) {
  donor <- as.character(donor)
  ug <- unique(donor)
  set.seed(seed)
  G <- nrow(counts)
  t_mat <- matrix(NA_real_, nrow = length(ug), ncol = G,
                  dimnames = list(ug, rownames(counts)))
  by_group <- split(seq_len(ncol(counts)), donor)
  K <- as.integer(K)

  for (d in ug) {
    ci <- by_group[[d]]
    Nd <- length(ci)
    if (Nd < 4L) next
    h <- Nd %/% 2L
    ## 2K sparse indicator columns, built as flat index vectors (no rbind loop)
    ii <- integer(0L); jj <- integer(0L)
    for (k in seq_len(K)) {
      perm <- sample.int(Nd)          # local positions inside this group
      ii <- c(ii, perm[seq_len(h)], perm[(h + 1L):Nd])
      jj <- c(jj, rep.int(2L * k - 1L, h), rep.int(2L * k, Nd - h))
    }
    W <- Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(Nd, 2L * K))
    s <- as.matrix(counts[, ci, drop = FALSE] %*% W)
    acc <- numeric(G); nk <- 0L
    for (k in seq_len(K)) {
      a <- s[, 2L * k - 1L]; b <- s[, 2L * k]
      sa <- sum(a); sb <- sum(b)
      if (sa <= 0 || sb <= 0) next
      yA <- log2(a / sa * 1e6 + 1)
      yB <- log2(b / sb * 1e6 + 1)
      acc <- acc + (yA - yB)^2 / 4
      nk <- nk + 1L
    }
    if (nk > 0L) t_mat[d, ] <- acc / nk
  }
  list(t_dg = t_mat, usable = rowSums(!is.na(t_mat)) > 0)
}
