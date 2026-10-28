#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 06_figure2.R
# Figure 2, the prescription, in the two forms a design question takes.
#
#   A  what a design can detect as donors are added, at several values of the
#      technical share rho
#   B  how many donors a target fold change costs, at the depth actually
#      observed, as a median with its interquartile range
#
#   MDE(n, rho) = (t_{1-a/2}(n-1) + t_{1-b}(n-1)) * sigma_bio * sqrt((1+rho)/n)
#
# Every gene keeps its own biological SD and, where the question is about this
# dataset, its own rho. The script asserts that it reproduces the two numbers
# the manuscript quotes (a 1.63-fold median at eight donors, and two-fold at six
# donors, IQR four to eight) and stops if it does not, rather than drawing a
# curve that disagrees with the text.
# ---------------------------------------------------------------------------
RES <- Sys.getenv("EFFN_RES", "C:/rtmp/effnfig2")
OUT <- Sys.getenv("EFFN_FIG", "C:/rtmp/effnfig2")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pg <- utils::read.delim(gzfile(file.path(RES, "gse96583_per_gene.tsv.gz")),
                        stringsAsFactors = FALSE)
ok <- pg[pg$status == "ok", ]
cat("estimable genes:", nrow(ok), "\n")

NMAX <- 60L
nn <- 2:NMAX
qn <- stats::qt(0.975, nn - 1) + stats::qt(0.80, nn - 1)     # design df = n - 1

curve_med <- function(n, rho) {
  stats::median(qn[n - 1] * ok$sigma_bio * sqrt((1 + rho) / n))
}

## donors each gene needs to reach a target fold change, at its own rho.
## qn is passed in rather than read from the enclosing frame so that the same
## routine serves both the nominal and the FDR-controlled threshold: the only
## difference between the two questions is the t quantile.
needed <- function(target, qn = NULL) {
  if (is.null(qn)) qn <- stats::qt(0.975, nn - 1) + stats::qt(0.80, nn - 1)
  lg <- log2(target)
  hit <- vapply(seq_len(nrow(ok)), function(i) {
    m <- qn * ok$sigma_bio[i] * sqrt((1 + ok$rho_tech[i]) / nn)
    j <- which(m <= lg)
    if (length(j)) nn[min(j)] else NA_real_
  }, 0.0)
  hit
}

## The FDR-controlled threshold uses the same arithmetic with an effective
## per-test level, so the two curves in panel B differ only in that level.
ALPHA_EFF <- 0.05 * 0.2 * 0.80 / ((1 - 0.2) * (1 - 0.05))
qn_fdr <- stats::qt(1 - ALPHA_EFF / 2, nn - 1) + stats::qt(0.80, nn - 1)
needed_fdr <- function(target) needed(target, qn_fdr)

focal <- stats::median(ok$rho_tech)
reported <- stats::median(ok$mde_fold)
anchor <- curve_med(8, ok$rho_tech)
cat(sprintf("median rho %.4f  median MDE at 8 donors %.6f  reported %.6f\n",
            focal, 2^anchor, reported))
stopifnot(abs(2^anchor - reported) < 1e-6)

n2 <- needed(2)
n15 <- needed(1.5)
q2 <- stats::quantile(n2, c(0.25, 0.5, 0.75), na.rm = TRUE)
cat(sprintf("two-fold: median %g (IQR %g to %g);  1.5-fold: median %g\n",
            q2[2], q2[1], q2[3], stats::median(n15, na.rm = TRUE)))
stopifnot(q2[2] == 6, q2[1] == 4, q2[3] == 8,
          stats::median(n15, na.rm = TRUE) == 11)

## the same two questions under 5% FDR control: the manuscript quotes eight
## donors for two-fold (IQR five to twelve) and sixteen for 1.5-fold, so the
## figure refuses to draw a curve that says anything else
n2f <- needed_fdr(2)
n15f <- needed_fdr(1.5)
q2f <- stats::quantile(n2f, c(0.25, 0.5, 0.75), na.rm = TRUE)
cat(sprintf("FDR: two-fold median %g (IQR %g to %g);  1.5-fold median %g\n",
            q2f[2], q2f[1], q2f[3], stats::median(n15f, na.rm = TRUE)))
stopifnot(q2f[2] == 8, q2f[1] == 5, q2f[3] == 12,
          stats::median(n15f, na.rm = TRUE) == 16)

## ---- panel A curves ------------------------------------------------------
ns <- 3:40
## rho = 1 is left out on purpose: the focal curve sits at 1.05 and the two
## would be drawn on top of one another
RHOS <- c(0.25, 0.5, 2, 4)
curves <- lapply(RHOS, function(r) vapply(ns, curve_med, 0.0, rho = r))
obs_curve <- vapply(ns, curve_med, 0.0, rho = focal)

## ---- panel B curve -------------------------------------------------------
tg <- seq(1.4, 3.0, by = 0.05)
bq <- t(vapply(tg, function(t) {
  v <- needed(t)
  stats::quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
}, c(0.0, 0.0, 0.0)))
bqf <- t(vapply(tg, function(t) {
  v <- needed_fdr(t)
  stats::quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE)
}, c(0.0, 0.0, 0.0)))

BLUE <- "#1f5fa9"; RED <- "#b8322a"; GREY <- "#6f6f6f"; DARK <- "#2f2f2f"
ramp <- c("#a9b8c9", "#8b9db3", "#63799a", "#3d5a80")

tiff(file.path(OUT, "Figure2_v200.tif"), width = 178, height = 91, units = "mm",
     res = 1200, compression = "lzw", type = "cairo")
par(mfrow = c(1, 2), mar = c(3.8, 4.1, 1.5, 0.9), mgp = c(2.4, 0.55, 0), tcl = -0.28,
    cex.axis = 0.85, cex.lab = 0.95, las = 1, col.axis = "#333333", fg = "#333333")

## ---- A ------------------------------------------------------------------
plot(NA, log = "xy", xlim = c(3, 40), ylim = c(1.1, 13),
     xlab = "Donors per group", ylab = "Median MDE (fold change)",
     xaxt = "n", yaxt = "n", bty = "l")
axis(1, at = c(3, 4, 6, 8, 12, 16, 24, 32, 40))
axis(2, at = c(1.2, 1.5, 2, 3, 5, 8, 12), labels = c(1.2, 1.5, 2, 3, 5, 8, 12))
for (i in seq_along(RHOS)) lines(ns, 2^curves[[i]], col = ramp[i], lwd = 1.5)
lines(ns, 2^obs_curve, col = BLUE, lwd = 2.6)
legend("topright", legend = c(sprintf("\u03c1 = %g", RHOS),
                              sprintf("\u03c1 = %.2f, this dataset", focal)),
       col = c(ramp, BLUE), lwd = c(rep(1.5, length(RHOS)), 2.6),
       bty = "n", cex = 0.7, seg.len = 1.5, inset = 0.015, y.intersp = 0.95)
abline(h = 2, col = GREY, lty = 2, lwd = 1.0)
text(40, 2, "two-fold", adj = c(1.06, -0.55), col = GREY, cex = 0.7)
abline(v = 8, col = RED, lty = 3, lwd = 1.0)
points(8, reported, pch = 21, bg = RED, col = RED, cex = 1.8, lwd = 1.3)
text(9.5, 2.95, sprintf("8 donors, %.2f-fold", reported),
     adj = c(0, 0), col = RED, cex = 0.74)
segments(8.15, reported * 1.06, 9.4, 2.8, col = RED, lwd = 0.6)
mtext("A", side = 3, line = 0.3, adj = -0.18, font = 2, cex = 1.0)

## ---- B ------------------------------------------------------------------
plot(NA, xlim = c(1.4, 3.0), ylim = c(0, 39), bty = "l",
     xlab = "Target fold change", ylab = "Donors required")
polygon(c(tg, rev(tg)), c(bqf[, 1], rev(bqf[, 3])), col = "#e8e5e0", border = NA)
polygon(c(tg, rev(tg)), c(bq[, 1], rev(bq[, 3])), col = "#dde5ef", border = NA)
lines(tg, bqf[, 2], col = DARK, lwd = 2.0, lty = 2)
lines(tg, bq[, 2], col = BLUE, lwd = 2.4)
abline(h = 8, col = RED, lty = 3, lwd = 1.0)
abline(v = 2, col = GREY, lty = 2, lwd = 1.0)
legend("topright", bty = "n", cex = 0.68, seg.len = 1.7, y.intersp = 1.18,
       legend = c(sprintf("two-fold, nominal: %g (IQR %g to %g)",
                          q2[2], q2[1], q2[3]),
                  sprintf("two-fold, 5%% FDR: %g (IQR %g to %g)",
                          q2f[2], q2f[1], q2f[3])),
       col = c(BLUE, DARK), lty = c(1, 2), lwd = c(2.4, 2.0))
text(2.98, 26, "red line: 8 donors in hand", adj = c(1, 0.5), col = RED, cex = 0.66)
mtext("B", side = 3, line = 0.3, adj = -0.18, font = 2, cex = 1.0)

dev.off()

cat("figure written:", file.path(OUT, "Figure2_v200.tif"), "\n")
cat(sprintf("at 8 donors, rho 0.25 to 4 spans %.2f to %.2f-fold\n",
            2^curves[[1]][6], 2^curves[[4]][6]))
