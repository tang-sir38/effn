# Internal helpers: pseudobulk normalisation and cell-level mean expression.

#' Pseudobulk a counts matrix to log2(CPM + 1)
#'
#' Sums raw counts within each level of `key` first, then normalises, so that
#' only one library-size factor enters the donor-level value.
#'
#' @param counts genes x cells matrix or Matrix.
#' @param key character vector of length ncol(counts); usually the donor, or
#'   the donor-by-group combination for paired designs.
#' @return levels(key) x genes matrix of log2(CPM + 1).
#' @keywords internal
#' @noRd
pb_norm <- function(counts, key) {
  key <- as.character(key)
  ug <- unique(key)
  out <- matrix(0, nrow = length(ug), ncol = nrow(counts),
                dimnames = list(ug, rownames(counts)))
  for (d in ug) {
    ci <- which(key == d)
    s <- Matrix::rowSums(counts[, ci, drop = FALSE])
    out[d, ] <- as.numeric(s)
  }
  L <- rowSums(out)
  if (any(L <= 0)) {
    stop("group with zero library size: ", paste(ug[L <= 0], collapse = ", "),
         call. = FALSE)
  }
  out <- log2(sweep(out, 1, L, "/") * 1e6 + 1)
  out
}

#' Cell-level mean expression in CPM (not logged)
#'
#' Computed as X %*% (1 / library size) / n_cells, i.e. a single sparse
#' matrix-vector product; no dense copy of the matrix is made.
#'
#' @param counts genes x cells matrix or Matrix.
#' @return numeric vector, one value per gene.
#' @keywords internal
#' @noRd
cell_mean_cpm <- function(counts) {
  L <- Matrix::colSums(counts)
  if (any(L <= 0)) stop("cell with zero library size", call. = FALSE)
  w <- 1 / L
  mu <- as.numeric(counts %*% w) * 1e6 / ncol(counts)
  names(mu) <- rownames(counts)
  mu
}

#' Column variances without the per-column R call overhead
#'
#' `apply(m, 2, var)` costs one interpreted call per gene, which dominates the
#' running time once a bootstrap loop is added. This is the same estimator,
#' vectorised.
#'
#' @param m numeric matrix, samples in rows.
#' @return numeric vector of column variances (denominator n - 1).
#' @keywords internal
#' @noRd
col_var <- function(m) {
  m <- as.matrix(m)
  n <- nrow(m)
  if (n < 2L) return(rep(NA_real_, ncol(m)))
  mu <- colMeans(m)
  (colSums(m * m) - n * mu * mu) / (n - 1)
}

#' Pooled within-group variance of donor-level values
#'
#' @param y donors x genes matrix.
#' @param group vector of length nrow(y) with exactly two distinct values.
#' @return list with `s2`, `df`, `n`.
#' @keywords internal
#' @noRd
sigma_pooled <- function(y, group) {
  group <- as.character(group)
  gs <- sort(unique(group))
  idx <- lapply(gs, function(k) which(group == k))
  n_k <- vapply(idx, length, integer(1))
  v_k <- lapply(idx, function(ii) col_var(y[ii, , drop = FALSE]))
  df <- sum(n_k - 1L)
  s2 <- Reduce(`+`, Map(function(v, n) v * (n - 1), v_k, n_k)) / df
  list(s2 = as.numeric(s2), df = df, n = stats::setNames(n_k, gs))
}

#' Design degrees of freedom for a donor-level contrast
#'
#' The t quantile belongs to the design, not to the variance ratio. A paired
#' design on `n` donors is a one-sample test on `n` within-donor differences and
#' carries `n - 1` degrees of freedom; an unpaired design on `n1` and `n2`
#' donors carries `n1 + n2 - 2`. The effective donor counts set the standard
#' error of the contrast, and they must not be used here: doing so shrinks the
#' degrees of freedom by the factor `1 + rho`, which inflates the t quantile
#' and therefore the MDE, exactly in the regime the package exists to describe.
#'
#' @keywords internal
#' @noRd
design_df <- function(n1, n2 = NULL, paired = FALSE) {
  if (isTRUE(paired)) return(rep_len(n1, length(n1)) - 1)
  rep_len(n1, length(n1)) + rep_len(n2, length(n1)) - 2
}

#' Minimum detectable effect from t quantiles
#'
#' Paired design: one-sample form on the within-donor contrast, with `df = n - 1`.
#' Unpaired design: two-sample form with equal allocation, with `df = n1 + n2 - 2`.
#'
#' `df` is required. It is the design degrees of freedom from [design_df()], not
#' the effective degrees of freedom `n_eff - 1`.
#'
#' @keywords internal
#' @noRd
mde_t <- function(sigma_bio, neff1, neff2, alpha = 0.05, power = 0.80,
                  paired = FALSE, df) {
  if (missing(df) || is.null(df)) {
    stop("mde_t() needs the design degrees of freedom (design_df()); the ",
         "effective donor counts enter the standard error, not the t quantile.",
         call. = FALSE)
  }
  df <- rep_len(df, length(neff1))
  sb <- rep_len(sigma_bio, length(df))
  out <- rep(NA_real_, length(df))
  ok <- is.finite(df) & df > 0 & is.finite(sb) & sb > 0 & is.finite(neff1)
  if (!any(ok)) return(out)
  tt <- stats::qt(1 - alpha / 2, df[ok]) + stats::qt(power, df[ok])
  if (paired) {
    out[ok] <- tt * sb[ok] * sqrt(1 / neff1[ok])
  } else {
    out[ok] <- tt * sb[ok] * sqrt(1 / neff1[ok] + 1 / neff2[ok])
  }
  out
}
