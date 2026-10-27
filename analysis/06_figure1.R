# Figure 1: three panels, all numbers read from 03_RESULTS (no hand-typed values).
#   A  per-gene effective donor count against the nominal eight
#   B  per-gene minimum detectable effect, in fold-change units
#   C  the depth verdict of the protocol: technical share against cells per donor
# Output: 04_FIGURES/Figure1.tif, 180 mm wide, 1200 dpi, LZW.
args <- commandArgs(trailingOnly = TRUE)
# The project path contains non-ASCII characters, and R decodes --file= in the
# native locale, so the derived path arrives mangled. Take the project root from
# the command line, or fall back to the working directory when the script is run
# from the project root.
root <- if (length(args) >= 1) args[1] else getwd()
stopifnot(dir.exists(file.path(root, "03_RESULTS")))
gse <- file.path(root, "03_RESULTS", "gse96583")
fig <- file.path(root, "04_FIGURES")
dir.create(fig, showWarnings = FALSE)

pg <- read.delim(gzfile(file.path(gse, "gse96583_per_gene.tsv.gz")), stringsAsFactors = FALSE)
dc <- read.delim(file.path(gse, "depth_curve.tsv"), stringsAsFactors = FALSE)
ok <- pg[pg$status == "ok", ]
stopifnot(nrow(ok) > 1000, !anyNA(ok$n_eff), !anyNA(ok$mde_log2))

neff <- ok$n_eff
med_neff <- median(neff)
med_fold <- median(2 ^ ok$mde_log2)
n_nominal <- 8

COL_BAR <- "#AFC3CC"
COL_ACC <- "#B03A2E"
COL_INK <- "#1A1A1A"

bars <- function(brk, cnt, ymax) {
  rect(brk[-length(brk)], 0, brk[-1], cnt, col = COL_BAR, border = NA)
  invisible(ymax)
}

tiff(file.path(fig, "Figure1.tif"), width = 180 / 25.4, height = 64 / 25.4,
     units = "in", res = 1200, compression = "lzw", type = "cairo")
# Cairo measures text against the device resolution, so the printed size is
# smaller than cex implies: at cex.axis = 0.80 a digit is 4.6 pt tall, below the
# 6 pt floor OUP applies to figure lettering. 0.95 puts it at about 7.5 pt.
par(mfrow = c(1, 3), mar = c(3.5, 4.1, 2.4, 0.8), mgp = c(2.0, 0.5, 0),
    tcl = -0.26, las = 1, cex.axis = 0.95, cex.lab = 1.0, ps = 12,
    col.axis = COL_INK, col.lab = COL_INK, family = "sans", lwd = 0.8)

## A --------------------------------------------------------------------------
br <- seq(1.9, 8.1, by = 0.1)
cnt <- hist(neff, breaks = br, plot = FALSE)$counts
ymax <- max(cnt) * 1.22
plot(NA, xlim = c(1.85, 8.15), ylim = c(0, ymax), axes = FALSE,
     xlab = expression(italic(n)[eff] * " per gene"), ylab = "genes")
bars(br, cnt, ymax)
axis(1, at = 2:8)
axis(2)
abline(v = n_nominal, lty = 3, lwd = 0.9)
abline(v = med_neff, lwd = 1.5, col = COL_ACC)
mtext("8 donors", side = 3, at = n_nominal, line = -0.35, cex = 0.79, adj = 1)
text(med_neff, ymax * 0.965, sprintf("%.2f", med_neff), cex = 0.95, col = COL_ACC,
     adj = c(0.5, 0.5), font = 2)
text(med_neff + 0.18, ymax * 0.90, "median", cex = 0.79, col = COL_ACC, adj = c(0, 0.5))
mtext("A", side = 3, adj = 0, line = 0.55, cex = 1.15, font = 2)

## B --------------------------------------------------------------------------
# fold change has a long right tail (21.9-fold), so bin on the log2 axis and
# label the ticks in fold units; genes above 6-fold are counted, not drawn.
hi <- log2(6)
n_above <- sum(ok$mde_log2 > hi)
br2 <- c(seq(0, hi, by = 0.05), hi + 1e-9)
cnt2 <- hist(pmin(ok$mde_log2, hi), breaks = br2, plot = FALSE)$counts
cnt2[length(cnt2)] <- 0   # the last bin is a pile-up of the clipped tail; count it in text
ymax2 <- max(cnt2) * 1.22
plot(NA, xlim = c(0, hi), ylim = c(0, ymax2), axes = FALSE,
     xlab = "minimum detectable effect (fold)", ylab = "genes")
bars(br2, cnt2, ymax2)
axis(1, at = log2(c(1, 1.5, 2, 3, 4, 6)), labels = c("1", "1.5", "2", "3", "4", "6"))
axis(2)
abline(v = log2(c(1.5, 2)), lty = 3, lwd = 0.8)
abline(v = log2(med_fold), lwd = 1.5, col = COL_ACC)
text(log2(med_fold), ymax2 * 0.965, sprintf("%.2f", med_fold), cex = 0.95, col = COL_ACC,
     adj = c(0.5, 0.5), font = 2)
text(log2(med_fold) - 0.09, ymax2 * 0.60, "median", cex = 0.79, col = COL_ACC, adj = c(1, 0.5))
if (n_above > 0) {
  # the tail is clipped at 6-fold, so say how much was clipped; the top right of
  # the panel is empty because the tail is short
  text(hi, ymax2 * 0.84, sprintf("%d genes above\n6-fold (not shown)", n_above),
       cex = 0.76, col = COL_INK, adj = c(1, 0))
}
mtext("B", side = 3, adj = 0, line = 0.55, cex = 1.15, font = 2)

## C --------------------------------------------------------------------------
x <- dc$cells_per_donor
y <- 100 * dc$core_depth_share
# xlim stops just past the last grid point: a padded panel leaves the axis line
# shorter than the reference line, which reads as a truncated plot
plot(NA, xlim = c(9, 1150), ylim = c(0, 100), log = "x", axes = FALSE,
     xlab = "cells per donor", ylab = "technical share of variance (%)")
# decade-style ticks; the last grid point (1,042) is named in the text and its
# label would collide with the one at 1,000
axis(1, at = c(10, 25, 100, 1000), labels = c("10", "25", "100", "1,000"),
     gap.axis = -1)
axis(2, at = seq(0, 100, 20))
abline(h = min(y), lty = 3, lwd = 0.9)
lines(x, y, lwd = 1.6, col = COL_ACC)
points(x, y, pch = 19, cex = 0.46, col = COL_ACC)
text(11.5, min(y) - 5.5, sprintf("lowest observed, %.0f%%", min(y)), cex = 0.76,
     adj = c(0, 0.5))
text(1042, y[length(y)] + 9, sprintf("%.0f%%", y[length(y)]), cex = 0.88,
     col = COL_ACC, adj = c(1, 0.5), font = 2)
text(x[1], y[1] + 13, sprintf("%.0f%%", y[1]), cex = 0.88, col = COL_ACC,
     adj = c(0, 0.5), font = 2)
mtext("C", side = 3, adj = 0, line = 0.55, cex = 1.15, font = 2)

invisible(dev.off())

cat("wrote", file.path(fig, "Figure1.tif"), "\n")
cat(sprintf("A: n_ok=%d median_n_eff=%.4f n_at_nominal=%d\n", nrow(ok), med_neff,
            sum(neff >= n_nominal - 1e-9)))
cat(sprintf("B: median_mde_fold=%.4f max=%.4f above_6fold=%d\n", med_fold,
            max(2 ^ ok$mde_log2), n_above))
cat(sprintf("C: %.1f%% at %d cells -> %.1f%% at %d cells; min %.1f%%; never below 20%%: %s\n",
            100 * dc$core_depth_share[1], x[1],
            100 * dc$core_depth_share[length(x)], x[length(x)], min(y),
            min(y) >= 20))
cat("verdict: technical share never below", sprintf("%.1f%%", min(y)), "\n")
