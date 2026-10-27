#' Plot an effn report
#'
#' Three panels: the distribution of the effective number of donors across
#' genes, the minimum detectable effect as a function of the number of donors
#' per group, and the loss ladder from cells to donors to effective donors.
#'
#' @param report An `effn_report` from [effn()].
#' @param which Character vector of panels to draw; any of `"neff"`, `"mde"`,
#'   `"loss"`.
#' @param target_mde Optional target effect (log2) drawn as a horizontal line in
#'   the MDE panel.
#' @param col Three colours, one per panel.
#' @param ... Passed to [graphics::par()].
#'
#' @return Invisibly, a list of the data frames used for each panel.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(2000, 5), nrow = 100,
#'                  dimnames = list(paste0("g", 1:100), NULL))
#' donor <- rep(paste0("d", 1:8), each = 25)
#' group <- rep(rep(c("ctrl", "case"), each = 4), each = 25)
#' rep <- effn(counts, group, donor, n_boot = 20, n_split = 5)
#' plot_effn(rep)
#'
#' @export
plot_effn <- function(report, which = c("neff", "mde", "loss"), target_mde = NULL,
                      col = c("#2B6CB0", "#C05621", "#2F855A"), ...) {
  stopifnot(inherits(report, "effn_report"))
  s <- report$summary
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  k <- length(which)
  graphics::par(mfrow = c(1, k), mar = c(4.2, 4.2, 2.2, 1), mgp = c(2.4, 0.7, 0), ...)
  out <- list()

  if ("neff" %in% which) {
    v <- report$per_gene$n_eff[report$per_gene$status == "ok"]
    v <- v[is.finite(v)]
    graphics::hist(v, breaks = 40, col = grDevices::adjustcolor(col[1], 0.35),
                   border = "white", main = "Effective donors per gene",
                   xlab = "n_eff", ylab = "genes")
    graphics::abline(v = s$n_donor, lty = 2, lwd = 1.4)
    graphics::abline(v = s$median_n_eff, lwd = 2, col = col[1])
    graphics::legend("topleft", bty = "n", cex = 0.8,
                     legend = c(sprintf("median %.2f", s$median_n_eff),
                                sprintf("donors %d", s$n_donor)),
                     col = c(col[1], "black"), lwd = c(2, 1.4), lty = c(1, 2))
    out$neff <- data.frame(n_eff = v)
  }

  if ("mde" %in% which) {
    n_cur <- min(s$n_donor_by_group)
    n_max <- max(ceiling(n_cur * 2.5), 12L)
    nn <- seq(3L, n_max)
    ne <- nn / (1 + s$median_rho_tech)
    m <- mde_t(s$median_sigma_bio, ne, ne,
               alpha = report$params$alpha, power = report$params$power,
               paired = isTRUE(report$params$paired),
               df = design_df(nn, nn, paired = isTRUE(report$params$paired)))
    graphics::plot(nn, m, type = "l", lwd = 2, col = col[2],
                   xlab = "donors per group", ylab = "MDE (log2)",
                   main = "Minimum detectable effect", ylim = c(0, max(m)))
    graphics::abline(v = n_cur, lty = 2)
    nce <- n_cur / (1 + s$median_rho_tech)
    graphics::points(n_cur, mde_t(s$median_sigma_bio, nce, nce,
                                  report$params$alpha, report$params$power,
                                  paired = isTRUE(report$params$paired),
                                  df = design_df(n_cur, n_cur,
                                                 paired = isTRUE(report$params$paired))),
                     pch = 19, col = col[2])
    if (!is.null(target_mde)) {
      graphics::abline(h = target_mde, lty = 3)
      n_req <- prescribe(report, target_mde = target_mde)$n_required
      if (!is.na(n_req)) graphics::abline(v = n_req, lty = 4, col = col[2])
    }
    out$mde <- data.frame(n_per_group = nn, mde_log2 = m)
  }

  if ("loss" %in% which) {
    v <- c(cells = s$nominal_n, donors = s$n_donor, n_eff = s$median_n_eff)
    graphics::barplot(log10(v), col = grDevices::adjustcolor(col[3], 0.45),
                      border = col[3], names.arg = c("cells", "donors", "n_eff"),
                      ylab = "log10(count)", main = "Loss ladder", las = 1)
    graphics::text(seq_along(v), log10(v) / 2, format(v, big.mark = ",", digits = 3),
                   cex = 0.75)
    out$loss <- data.frame(step = names(v), value = as.numeric(v))
  }

  invisible(out)
}
