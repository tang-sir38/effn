#' Effective sample size, minimum detectable effect and prescription
#'
#' Takes a single-cell count matrix with donor labels and returns, per gene,
#' the effective number of donors after the depth penalty, the minimum
#' detectable effect at the requested power, and how many additional donors a
#' target effect would require.
#'
#' @param object A genes x cells count matrix (`matrix`, `Matrix`,
#'   `data.frame`), a `SingleCellExperiment`, or a Seurat object (v4 or v5).
#' @param group Vector of length `ncol(object)` giving the group of each cell.
#'   Exactly two groups are allowed.
#' @param donor Optional vector of donor labels. When `NULL`, labels are looked
#'   up in sample metadata; if none are found the function stops rather than
#'   treating the matrix as a single replicate.
#' @param assay Assay name for `SingleCellExperiment` / Seurat input.
#' @param layer Layer (Seurat v5) or slot (Seurat v4) holding raw counts.
#' @param n_boot Number of donor-level bootstrap replicates.
#' @param n_tech Number of within-donor resamples used to measure the technical
#'   variance of the pseudobulk mean.
#' @param seed Integer seed; identical input and seed give identical output.
#' @param alpha Two-sided type I error.
#' @param power Target power.
#' @param min_donor Minimum donors per group; below this the gene is reported as
#'   unestimable rather than given a number.
#' @param tech `"bootstrap"` (default) resamples each donor's cells with
#'   replacement at the observed depth; `"splithalf"` uses the cheaper split-half
#'   estimator, whose 1/n divisor only holds on a linear scale, not on
#'   log2(CPM + 1).
#' @param paired `"auto"` (default) treats the design as paired when every
#'   donor contributes to both groups; `TRUE` / `FALSE` force the choice.
#' @param jack Use a leave-one-donor-out jackknife for the reported intervals of
#'   the median effective count and the median MDE (paired designs). The donor
#'   bootstrap duplicates donors and shrinks the between-donor variance, so its
#'   interval sits below the point estimate; the jackknife does not.
#' @param verbose Print stage timings to the message stream.
#'
#' @return An object of class `effn_report` with elements `per_gene`,
#'   `summary`, `boot`, `params` and `call`.
#'
#' @examples
#' set.seed(1)
#' G <- 60; d <- rep(paste0("d", 1:8), each = 20)
#' g <- rep(rep(c("ctrl", "case"), each = 4), each = 20)
#' mu <- exp(rnorm(G, 2, 0.6)); b <- matrix(rnorm(G * 8, 0, 0.5), G)
#' cnt <- sapply(seq_along(d), function(j)
#'   rnbinom(G, mu = mu * exp(b[, match(d[j], paste0("d", 1:8))]), size = 10))
#' rownames(cnt) <- paste0("g", 1:G)
#' r <- effn(cnt, g, d, n_boot = 20, n_tech = 20)
#' summary(r)
#'
#' @export
effn <- function(object, group, donor = NULL, assay = NULL, layer = "counts",
                 n_boot = 200L, n_tech = 40L,
                 tech = c("bootstrap", "splithalf"), seed = 20261021L,
                 alpha = 0.05, power = 0.80, min_donor = 3L, paired = "auto",
                 jack = TRUE, verbose = FALSE) {

  tech <- match.arg(tech)

  stamp <- function(...) if (isTRUE(verbose)) {
    message(format(Sys.time(), "[%H:%M:%S] "), ...)
  }

  counts <- as_counts(object, assay = assay, layer = layer)
  ## readMM() and many on-disk formats hand back triplet matrices; every step
  ## below indexes columns, so normalise to compressed sparse once (this also
  ## sums repeated (i, j) entries).
  if (inherits(counts, "TsparseMatrix")) counts <- methods::as(counts, "CsparseMatrix")
  if (!is.numeric(counts[1, 1]) && !inherits(counts, "Matrix")) {
    stop("counts must be numeric", call. = FALSE)
  }
  if (is.null(rownames(counts))) rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell", seq_len(ncol(counts)))

  ## "group" may be the name of a metadata column rather than a vector
  group <- resolve_column(object, group)

  d <- extract_donor(object, donor)
  donor <- d$donor
  if (length(donor) != ncol(counts)) {
    stop("length(donor) = ", length(donor), " but ncol(counts) = ", ncol(counts),
         call. = FALSE)
  }
  if (length(group) != ncol(counts)) {
    stop("length(group) = ", length(group), " but ncol(counts) = ", ncol(counts),
         call. = FALSE)
  }
  group <- as.character(group)
  gs <- sort(unique(group))
  if (length(gs) != 2L) {
    stop("exactly two groups are required, got ", length(gs), ": ",
         paste(gs, collapse = ", "), call. = FALSE)
  }

  tab <- table(donor, group)
  n_in <- rowSums(tab > 0)
  if (identical(paired, "auto")) {
    paired <- all(n_in == 2L) && length(unique(donor)) >= min_donor
  }
  if (isTRUE(paired)) {
    if (!all(n_in == 2L)) {
      stop("paired = TRUE requires every donor in both groups; donors with one ",
           "group: ", paste(utils::head(names(n_in)[n_in != 2L], 5), collapse = ", "),
           call. = FALSE)
    }
    n_k <- stats::setNames(rep(length(unique(donor)), 2L), gs)
  } else {
    if (any(n_in > 1L)) {
      stop("donor(s) appear in both groups; set paired = TRUE for a paired design. ",
           "Offending: ", paste(utils::head(names(n_in)[n_in > 1L], 5), collapse = ", "),
           call. = FALSE)
    }
    n_k <- vapply(gs, function(k) length(unique(donor[group == k])), integer(1))
  }
  if (any(n_k < min_donor)) {
    stop("each group needs at least ", min_donor, " donors; got ",
         paste(paste0(names(n_k), "=", n_k), collapse = ", "), call. = FALSE)
  }

  ## genes that are all zero are reported, not silently dropped
  tot <- as.numeric(Matrix::rowSums(counts))
  all_zero <- tot <= 0
  keep <- !all_zero
  counts_k <- counts[keep, , drop = FALSE]

  ## donor-by-group key: one pseudobulk value per (donor, group) cell
  key <- paste(donor, group, sep = "\r")
  y <- pb_norm(counts_k, key)
  stamp(nrow(counts_k), " genes x ", ncol(counts_k), " cells, ",
        nrow(y), " donor-group pseudobulks")
  kd <- sub("\r.*$", "", rownames(y))
  kg <- sub("^.*\r", "", rownames(y))

  if (isTRUE(paired)) {
    ## wide layout: one row per donor, columns ctrl / case
    ord <- order(kd, kg)
    y_w <- y[ord, , drop = FALSE]
    dd <- y_w[kg[ord] == gs[2], , drop = FALSE] - y_w[kg[ord] == gs[1], , drop = FALSE]
    rownames(dd) <- kd[ord][kg[ord] == gs[2]]
    s2_obs <- col_var(dd)
    df_obs <- nrow(dd) - 1L
  } else {
    sp <- sigma_pooled(y, kg)
    s2_obs <- sp$s2
    df_obs <- sp$df
  }

  tv <- if (identical(tech, "bootstrap")) {
    sigma_tech_boot(counts_k, key, B = n_tech, seed = seed)
  } else {
    sigma_tech_splithalf(counts_k, key, K = n_tech, seed = seed)
  }
  stamp("technical variance done (", n_tech, " ", tech, " resamples per donor-group)")
  t_mat <- tv$t_dg
  if (isTRUE(paired)) {
    tk <- sub("\r.*$", "", rownames(t_mat))
    tg <- sub("^.*\r", "", rownames(t_mat))
    t1 <- t_mat[tg == gs[1], , drop = FALSE]
    t2 <- t_mat[tg == gs[2], , drop = FALSE]
    rownames(t1) <- tk[tg == gs[1]]; rownames(t2) <- tk[tg == gs[2]]
    t2 <- t2[rownames(t1), , drop = FALSE]
    t_donor <- t1 + t2                       # variance of the difference
  } else {
    t_donor <- t_mat
  }
  t_use <- t_donor[rowSums(!is.na(t_donor)) > 0, , drop = FALSE]
  s2_tech <- if (nrow(t_use) > 0L) colMeans(t_use) else rep(NA_real_, length(s2_obs))

  mu <- cell_mean_cpm(counts_k)
  s2_bio <- s2_obs - s2_tech
  tech_dominated <- is.na(s2_bio) | s2_bio <= 0
  s2_bio_c <- ifelse(tech_dominated, 0, s2_bio)
  rho <- ifelse(s2_bio_c > 0, s2_tech / s2_bio_c, NA_real_)
  n1 <- n_k[[1]]; n2 <- n_k[[2]]
  neff1 <- n1 / (1 + rho); neff2 <- n2 / (1 + rho)
  if (any(neff1 > n1 + 1e-8, na.rm = TRUE) || any(neff2 > n2 + 1e-8, na.rm = TRUE)) {
    stop("internal error: n_eff exceeds n_donor (rho_tech < 0)", call. = FALSE)
  }

  mde <- rep(NA_real_, length(s2_obs))
  ## Degrees of freedom come from the design. A resolution floor on the
  ## effective count is a separate decision: below two effective donors the
  ## contrast rests on a single replicate, so we withhold the number rather
  ## than report a large but finite MDE.
  df_design <- design_df(n1, n2, paired = isTRUE(paired))
  res_ok <- if (isTRUE(paired)) (neff1 - 1 > 0) else (neff1 + neff2 - 2 > 0)
  mdx <- which(!is.na(neff1) & !is.na(neff2) & neff1 > 0 & res_ok)
  if (length(mdx) > 0L) {
    mde[mdx] <- mde_t(sqrt(s2_bio_c)[mdx], neff1[mdx], neff2[mdx],
                      alpha = alpha, power = power, paired = isTRUE(paired),
                      df = df_design)
  }

  status_k <- rep("ok", length(s2_obs))
  status_k[mu < 0.1] <- "low_expression"
  status_k[tech_dominated] <- "tech_dominated"
  ## Resolution floor: with fewer than two effective donors the t reference
  ## distribution has under one degree of freedom in a paired design, and the
  ## quantile ratio explodes (observed ratio > 1e4 on GSE96583). Such a gene has
  ## less information than a two-donor comparison, so it is reported as
  ## low_resolution rather than given a number that is really an artefact of a
  ## vanishing denominator.
  low_res <- !tech_dominated & is.finite(neff1) & neff1 < 2
  status_k[low_res] <- "low_resolution"
  mde[low_res] <- NA_real_
  status_k[status_k == "ok" & is.na(mde)] <- "uninformative"

  pg <- data.frame(
    gene = rownames(counts),
    mu_cpm = NA_real_, mean_log2cpm = NA_real_,
    sigma_obs = NA_real_, sigma_tech = NA_real_, sigma_bio = NA_real_,
    rho_tech = NA_real_, n_eff = NA_real_,
    n_eff_group1 = NA_real_, n_eff_group2 = NA_real_,
    mde_log2 = NA_real_, mde_fold = NA_real_,
    status = "all_zero", stringsAsFactors = FALSE
  )
  pg$mu_cpm[keep] <- mu
  pg$mean_log2cpm[keep] <- colMeans(y)
  pg$sigma_obs[keep] <- sqrt(s2_obs)
  pg$sigma_tech[keep] <- sqrt(s2_tech)
  pg$sigma_bio[keep] <- sqrt(s2_bio_c)
  pg$rho_tech[keep] <- rho
  pg$n_eff[keep] <- neff1
  pg$n_eff_group1[keep] <- neff1
  pg$n_eff_group2[keep] <- neff2
  pg$mde_log2[keep] <- mde
  pg$mde_fold[keep] <- 2^mde
  pg$status[keep] <- status_k
  per_gene <- pg

  ## ---- bootstrap over donors -------------------------------------------
  ok <- which(status_k == "ok")
  boot_neff <- boot_mde <- boot_sig <- boot_rho <- rep(NA_real_, n_boot)
  if (length(ok) > 0L && n_boot > 0L) {
    set.seed(seed)
    if (isTRUE(paired)) {
      donors_all <- rownames(dd)
      for (b in seq_len(n_boot)) {
        ii <- sample(donors_all, length(donors_all), replace = TRUE)
        sb <- col_var(dd[ii, , drop = FALSE])
        tt <- colMeans(t_donor[ii, , drop = FALSE], na.rm = TRUE)
        bio <- sb - tt
        r <- ifelse(bio > 0, tt / bio, NA_real_)
        kk <- ok[is.finite(r[ok]) & r[ok] >= 0]
        boot_neff[b] <- stats::median(n1 / (1 + r[kk]), na.rm = TRUE)
        if (length(kk) > 0L) {
          boot_mde[b] <- stats::median(
            mde_t(sqrt(pmax(bio[kk], 0)), n1 / (1 + r[kk]), n2 / (1 + r[kk]),
                  alpha = alpha, power = power, paired = TRUE, df = df_design),
            na.rm = TRUE)
        }
        boot_sig[b] <- stats::median(sqrt(pmax(bio[kk], 0)), na.rm = TRUE)
        boot_rho[b] <- stats::median(r[kk], na.rm = TRUE)
      }
    } else {
      d1 <- unique(donor[group == gs[1]]); d2 <- unique(donor[group == gs[2]])
      for (b in seq_len(n_boot)) {
        i1 <- sample(d1, length(d1), replace = TRUE)
        i2 <- sample(d2, length(d2), replace = TRUE)
        rows <- c(paste0(i1, "\r", gs[1]), paste0(i2, "\r", gs[2]))
        yy <- y[rows, , drop = FALSE]
        gg <- c(rep(gs[1], length(i1)), rep(gs[2], length(i2)))
        sb <- sigma_pooled(yy, gg)$s2
        tt <- colMeans(t_donor[rows, , drop = FALSE], na.rm = TRUE)
        bio <- sb - tt
        r <- ifelse(bio > 0, tt / bio, NA_real_)
        kk <- ok[is.finite(r[ok]) & r[ok] >= 0]
        boot_neff[b] <- stats::median(n1 / (1 + r[kk]), na.rm = TRUE)
        if (length(kk) > 0L) {
          boot_mde[b] <- stats::median(
            mde_t(sqrt(pmax(bio[kk], 0)), n1 / (1 + r[kk]), n2 / (1 + r[kk]),
                  alpha = alpha, power = power, df = df_design), na.rm = TRUE)
        }
        boot_sig[b] <- stats::median(sqrt(pmax(bio[kk], 0)), na.rm = TRUE)
        boot_rho[b] <- stats::median(r[kk], na.rm = TRUE)
      }
    }
  }
  ci <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 10L) return(c(NA_real_, NA_real_))
    stats::quantile(x, c(0.025, 0.975), names = FALSE)
  }
  ci_neff <- ci(boot_neff); ci_mde <- ci(boot_mde)
  stamp("bootstrap done (", n_boot, " donor resamples)")

  ## Reported intervals come from a leave-one-donor-out jackknife when the
  ## design allows it. The donor bootstrap resamples with replacement, so
  ## duplicated donors enter the between-donor variance and shrink it; the
  ## resulting interval sits systematically below the point estimate and can
  ## fail to cover it. Deleting one donor at a time does not duplicate anyone.
  pt_neff <- stats::median((n1 / (1 + rho))[ok], na.rm = TRUE)
  pt_mde <- stats::median(mde[ok], na.rm = TRUE)
  ci_method <- "donor bootstrap"
  jk_store <- NULL
  if (isTRUE(jack) && isTRUE(paired)) {
    nd <- nrow(dd)
    jn <- jm <- rep(NA_real_, nd)
    for (i in seq_len(nd)) {
      keep <- setdiff(seq_len(nd), i)
      s2i <- col_var(dd[keep, , drop = FALSE])
      ti <- colMeans(t_donor[keep, , drop = FALSE], na.rm = TRUE)
      bi <- s2i - ti
      ri <- ifelse(bi > 0, ti / bi, NA_real_)
      kk <- ok[is.finite(ri[ok]) & ri[ok] >= 0]
      if (!length(kk)) next
      jn[i] <- stats::median((n1 / (1 + ri))[kk], na.rm = TRUE)
      jm[i] <- stats::median(
        mde_t(sqrt(pmax(bi[kk], 0)), n1 / (1 + ri[kk]), n2 / (1 + ri[kk]),
              alpha = alpha, power = power, paired = TRUE, df = df_design),
        na.rm = TRUE)
    }
    jk_ci <- function(point, rp) {
      v <- rp[is.finite(rp)]
      if (length(v) < 3L) return(c(NA_real_, NA_real_))
      se <- sqrt((length(v) - 1) / length(v) * sum((v - mean(v))^2))
      hw <- stats::qt(1 - alpha / 2, length(v) - 1) * se
      c(point - hw, point + hw)
    }
    ci_neff <- jk_ci(pt_neff, jn)
    ci_mde <- jk_ci(pt_mde, jm)
    ci_method <- "jackknife (delete one donor)"
    jk_store <- list(n_eff = jn, mde = jm)
    stamp("jackknife done (", nd, " delete-one passes)")
  }

  n_nom <- ncol(counts)
  n_donor <- length(unique(donor))
  smed <- function(z) suppressWarnings(stats::median(z, na.rm = TRUE))
  rho_med <- smed(rho[ok])
  sig_med <- smed(sqrt(s2_bio_c)[ok])

  sm <- list(
    design = if (isTRUE(paired)) "paired" else "unpaired",
    tech_estimator = tech,
    nominal_n = n_nom,
    nominal_unit = "cells",
    n_donor = n_donor,
    n_donor_by_group = stats::setNames(as.integer(n_k), gs),
    n_genes_tested = nrow(per_gene),
    median_n_eff = pt_neff,
    median_n_eff_ci = ci_neff,
    median_mde_log2 = pt_mde,
    median_mde_fold = 2^pt_mde,
    median_mde_log2_ci = ci_mde,
    ci_method = ci_method,
    median_sigma_bio = sig_med,
    median_rho_tech = rho_med,
    loss_donor = n_nom / n_donor,
    loss_depth = 1 + rho_med,
    loss_depth_share = rho_med / (1 + rho_med),
    n_unestimable = c(
      all_zero = sum(all_zero),
      tech_dominated = sum(per_gene$status == "tech_dominated"),
      low_resolution = sum(per_gene$status == "low_resolution"),
      low_expression = sum(per_gene$status == "low_expression"),
      uninformative = sum(per_gene$status == "uninformative")
    ),
    n_ok = length(ok),
    groups = gs,
    donor_source = d$source,
    heuristic = d$heuristic
  )
  nreq2 <- prescribe_n_vec(1, sqrt(s2_bio_c)[ok], rho[ok], alpha = alpha,
                           power = power, paired = isTRUE(paired))
  nreq15 <- prescribe_n_vec(log2(1.5), sqrt(s2_bio_c)[ok], rho[ok],
                            alpha = alpha, power = power, paired = isTRUE(paired))
  n_cur <- min(n_k)
  sm$n_required_2fold <- prescribe_median(nreq2)
  sm$n_required_2fold_iqr <- if (all(is.na(nreq2))) c(NA_real_, NA_real_) else
    as.numeric(stats::quantile(nreq2, c(0.25, 0.75), na.rm = TRUE, names = FALSE))
  sm$n_required_1.5fold <- prescribe_median(nreq15)
  sm$n_required_n_unreachable <- sum(is.na(nreq2))
  sm$shortfall_2fold <- if (is.na(sm$n_required_2fold)) NA_integer_ else
    max(0L, as.integer(ceiling(sm$n_required_2fold)) - as.integer(n_cur))
  sm$prescription <- build_sentence(1, power, sm$n_required_2fold, n_cur,
                                    sm$shortfall_2fold, isTRUE(paired))
  sm$prescription_1.5fold <- build_sentence(log2(1.5), power, sm$n_required_1.5fold,
                                            n_cur,
                                            if (is.na(sm$n_required_1.5fold)) NA_integer_ else
                                              max(0L, as.integer(ceiling(sm$n_required_1.5fold)) - as.integer(n_cur)),
                                            isTRUE(paired))
  boot <- list(n_eff = boot_neff, mde = boot_mde, sigma_bio = boot_sig, rho = boot_rho)

  structure(list(per_gene = per_gene, summary = sm, boot = boot, jack = jk_store,
                 params = list(alpha = alpha, power = power, n_boot = n_boot,
                               n_tech = n_tech, tech = tech, seed = seed,
                               paired = isTRUE(paired), jack = isTRUE(jack),
                               groups = gs, min_donor = min_donor),
                 call = match.call()),
            class = "effn_report")
}
