#' Prescribe the number of donors needed for a target effect
#'
#' Inverts the minimum detectable effect formula: given the biological SD and
#' the depth penalty of a gene, returns how many donors per group are required
#' to detect `target_mde` at the requested power.
#'
#' With `gene = NULL` the returned number is the median over estimable genes of
#' the donor count that each gene would need. That is the same "typical gene"
#' summary that [effn()] reports for the MDE; inverting at the median SD instead
#' would answer a different question and disagree with the reported MDE whenever
#' the SD distribution is skewed.
#'
#' @param report An `effn_report` from [effn()].
#' @param target_mde Target effect on the log2 scale (e.g. `1` is a two-fold
#'   change, `log2(1.5)` a 1.5-fold change).
#' @param gene Optional gene name; when `NULL` the median across estimable genes
#'   is used.
#' @param alpha Two-sided type I error.
#' @param power Target power.
#' @param nmax Largest donor count considered; genes that would need more are
#'   reported as unreachable instead of being given an extrapolated number.
#'
#' @return A list with `n_required`, `n_required_ci`, `n_current`, `shortfall`,
#'   `target_mde`, `target_fold`, `paired`, `n_unreachable` and `sentence`.
#'
#' @examples
#' set.seed(1)
#' G <- 60; d <- rep(paste0("d", 1:8), each = 20)
#' g <- rep(rep(c("ctrl", "case"), each = 4), each = 20)
#' mu <- exp(rnorm(G, 2, 0.6)); b <- matrix(rnorm(G * 8, 0, 0.5), G)
#' cnt <- sapply(seq_along(d), function(j)
#'   rnbinom(G, mu = mu * exp(b[, match(d[j], paste0("d", 1:8))]), size = 10))
#' rownames(cnt) <- paste0("g", 1:G)
#' r <- effn(cnt, g, d, n_boot = 20, n_split = 5)
#' prescribe(r, target_mde = 1)
#'
#' @export
prescribe <- function(report, target_mde = 1, gene = NULL,
                      alpha = report$params$alpha, power = report$params$power,
                      nmax = 10000L) {
  stopifnot(inherits(report, "effn_report"))
  sm <- report$summary
  paired <- isTRUE(report$params$paired)
  n_cur <- min(sm$n_donor_by_group)

  if (!is.null(gene)) {
    row <- report$per_gene[report$per_gene$gene == gene, , drop = FALSE]
    if (nrow(row) == 0L) stop("gene not found: ", gene, call. = FALSE)
    if (!identical(row$status, "ok")) {
      return(list(n_required = NA_real_, n_required_ci = c(NA_real_, NA_real_),
                  n_current = as.integer(n_cur), shortfall = NA_integer_,
                  target_mde = target_mde, target_fold = 2^target_mde,
                  paired = paired, n_unreachable = NA_integer_,
                  sentence = paste0("Gene ", gene, " is ", row$status,
                                    "; no prescription is given.")))
    }
    n_req <- as.numeric(prescribe_n_vec(target_mde, row$sigma_bio, row$rho_tech,
                                        alpha = alpha, power = power,
                                        paired = paired, nmax = nmax))
    n_unreach <- if (is.na(n_req)) 1L else 0L
  } else {
    okg <- report$per_gene$status == "ok"
    n_all <- prescribe_n_vec(target_mde, report$per_gene$sigma_bio[okg],
                             report$per_gene$rho_tech[okg],
                             alpha = alpha, power = power, paired = paired,
                             nmax = nmax)
    n_req <- prescribe_median(n_all)
    n_unreach <- sum(is.na(n_all))
  }

  ## Interval from the donor-level bootstrap, applied to the resampled median
  ## SD and depth penalty. It is indicative for the median-over-genes number
  ## above, not a quantile of it.
  n_ci <- c(NA_real_, NA_real_)
  if (!is.null(report$boot) && length(report$boot$sigma_bio) > 10L) {
    bs <- report$boot$sigma_bio
    br <- report$boot$rho
    okb <- is.finite(bs) & is.finite(br) & bs > 0 & is.finite(br) & br >= 0
    if (sum(okb) >= 10L) {
      nr <- prescribe_n_vec(target_mde, bs[okb], br[okb], alpha = alpha,
                            power = power, paired = paired, nmax = nmax)
      nr <- nr[is.finite(nr)]
      if (length(nr) >= 10L) n_ci <- stats::quantile(nr, c(0.025, 0.975), names = FALSE)
    }
  }

  short <- if (is.na(n_req)) NA_integer_ else
    max(0L, as.integer(ceiling(n_req)) - as.integer(n_cur))
  sent <- build_sentence(target_mde, power, n_req, n_cur, short, paired)

  list(n_required = n_req, n_required_ci = n_ci, n_current = as.integer(n_cur),
       shortfall = short, target_mde = target_mde, target_fold = 2^target_mde,
       paired = paired, n_unreachable = n_unreach, sentence = sent)
}

#' Solve for the number of donors per group
#'
#' The MDE falls strictly as `n` grows, so the smallest `n` that reaches
#' `target_mde` is found by bisection. A direct fixed-point iteration on
#' `n` does not converge: the t quantiles depend on `n` through the degrees of
#' freedom, and the map can oscillate between two values indefinitely.
#'
#' @return integer vector, `NA` where the target is out of reach within `nmax`.
#' @keywords internal
#' @noRd
prescribe_n_vec <- function(target_mde, sigma_bio, rho, alpha = 0.05,
                            power = 0.80, nmax = 10000L, paired = FALSE,
                            n_min = 3L) {
  sigma_bio <- as.numeric(sigma_bio)
  g <- length(sigma_bio)
  rho <- rep_len(as.numeric(rho), g)
  out <- rep(NA_integer_, g)
  if (g == 0L) return(out)
  usable <- is.finite(sigma_bio) & sigma_bio > 0 & is.finite(rho) & rho >= 0 &
    is.finite(target_mde) & target_mde > 0
  if (!any(usable)) return(out)

  mde_at <- function(n, idx) {
    r <- rho[idx]
    sb <- sigma_bio[idx]
    ne <- n / (1 + r)
    ## design degrees of freedom for the candidate design: a paired study on n
    ## donors has n - 1, an unpaired one with n per group has 2n - 2
    df <- design_df(n, if (paired) NULL else n, paired = paired)
    tt <- suppressWarnings(stats::qt(1 - alpha / 2, df) + stats::qt(power, df))
    bad <- !is.finite(tt) | !is.finite(df) | df <= 0
    m <- if (paired) tt * sb * sqrt(1 / ne) else tt * sb * sqrt(2 / ne)
    m[bad] <- Inf
    as.numeric(m)
  }

  hi0 <- as.integer(nmax)
  reach <- mde_at(rep(hi0, g), seq_len(g)) <= target_mde
  usable <- usable & reach
  if (!any(usable)) return(out)

  lo <- rep(as.integer(n_min), g)
  hi <- rep(hi0, g)
  lo[!usable] <- NA_integer_
  hi[!usable] <- NA_integer_
  for (it in seq_len(64L)) {
    act <- which(!is.na(lo) & lo < hi)
    if (length(act) == 0L) break
    mid <- (lo[act] + hi[act]) %/% 2L
    go <- mde_at(mid, act) <= target_mde
    go[is.na(go)] <- FALSE
    hi[act[go]] <- mid[go]
    lo[act[!go]] <- mid[!go] + 1L
  }
  out[usable] <- lo[usable]
  out
}

#' @keywords internal
#' @noRd
prescribe_n <- function(target_mde, sigma_bio, rho, n_start = 3L,
                        alpha = 0.05, power = 0.80, nmax = 10000L,
                        paired = FALSE) {
  as.numeric(prescribe_n_vec(target_mde, sigma_bio, rho, alpha = alpha,
                             power = power, nmax = nmax, paired = paired,
                             n_min = n_start))
}

#' Median donor requirement, refusing to answer when most genes are unreachable
#' @keywords internal
#' @noRd
prescribe_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  as.numeric(stats::median(x))
}

#' One-sentence prescription
#' @keywords internal
#' @noRd
build_sentence <- function(target_mde, power, n_req, n_cur, short, paired) {
  unit <- if (paired) "donors (paired design)" else "donors per group"
  if (length(n_req) == 0L || !is.finite(n_req)) {
    return(sprintf(paste0("A %.1f-fold change is not reachable within 10,000 %s at the ",
                          "observed biological SD."), 2^target_mde, unit))
  }
  n_req <- as.integer(ceiling(n_req))
  if (short == 0L) {
    return(sprintf(paste0("Detecting a %.1f-fold change at %d%% power requires %d %s; ",
                          "this dataset has %d, so the design already supports it."),
                   2^target_mde, round(100 * power), as.integer(n_req), unit,
                   as.integer(n_cur)))
  }
  sprintf(paste0("Detecting a %.1f-fold change at %d%% power requires %d %s; ",
                 "this dataset has %d, %d short."),
          2^target_mde, round(100 * power), as.integer(n_req), unit,
          as.integer(n_cur), as.integer(short))
}

#' Sentence built from the median SD and depth penalty
#' @keywords internal
#' @noRd
prescribe_sentence <- function(sm, target_mde = 1, paired = FALSE,
                               alpha = 0.05, power = 0.80, nmax = 10000L) {
  n_req <- prescribe_n_vec(target_mde, sm$median_sigma_bio, sm$median_rho_tech,
                           alpha = alpha, power = power, paired = paired,
                           nmax = nmax)
  n_cur <- min(sm$n_donor_by_group)
  short <- if (is.na(n_req)) NA_integer_ else
    max(0L, as.integer(n_req) - as.integer(n_cur))
  build_sentence(target_mde, power, as.numeric(n_req), n_cur, short, paired)
}
