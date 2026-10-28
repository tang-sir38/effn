#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 06_figure1_v120.R
# 图 1 四面板。与 v1.1.2 的差别：
#   A  不变（n_eff 直方图），但标注中位数
#   B  加一条 FDR 校正后的竖线（1.93），原来的 1.63 保留为实线
#   C  横轴改为"每条件每供体细胞数"，纵轴改为 MDE（原来只有技术占比）
#   D  新增：同一批基因的 n_eff 随深度变化，用来说明 n_eff 不等于精度
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
## 本机 R 在中文路径下 gzfile() 会报 "unable to translate ... to native
## encoding"（getwd() 返回 <U+5510> 这类转义），normalizePath 又会展开成 8.3
## 短名。所以输入输出都走 ASCII 目录，由 EFFN_RES / EFFN_FIG 指定，
## 生成后再由外层脚本拷回项目目录。
RES <- Sys.getenv("EFFN_RES", "C:/rtmp/effnfig")
OUT <- Sys.getenv("EFFN_FIG", "C:/rtmp/effnfig")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

pg <- utils::read.delim(gzfile(file.path(RES, "gse96583_per_gene.tsv.gz")),
                        stringsAsFactors = FALSE)
ok <- pg[pg$status == "ok", ]
dc <- utils::read.delim(file.path(RES, "depth_curve_per_condition.tsv"))
dc$cells_per_donor_per_condition[is.na(dc$cells_per_donor_per_condition)] <- 1679

med_neff <- stats::median(ok$n_eff)
med_mde <- stats::median(ok$mde_fold)
mde_fdr <- 1.928                      # 03_RESULTS/gse96583/fdr_mde_summary.tsv, pi1 = 0.2

BLUE <- "#1f5fa9"; RED <- "#b8322a"; GREY <- "#6f6f6f"

pan <- function() {
  par(mar = c(3.6, 3.8, 1.6, 0.8), mgp = c(2.2, 0.5, 0), tcl = -0.25,
      cex.axis = 0.86, cex.lab = 0.95, las = 1, col.axis = "#333333",
      fg = "#333333", lend = "butt")
}

tiff(file.path(OUT, "Figure1_v200.tif"), width = 178, height = 148, units = "mm",
     res = 1200, compression = "lzw", type = "cairo")
par(mfrow = c(2, 2), oma = c(0, 0, 0, 0))

## ---- A: 有效供体数 ------------------------------------------------------
pan()
h <- hist(ok$n_eff, breaks = seq(0.5, 8.5, by = 0.25), plot = FALSE)
plot(h, col = "#cfdcec", border = NA, main = "", xlab = "Effective donors per gene",
     ylab = "Genes", xlim = c(0, 8.4), ylim = c(0, max(h$counts) * 1.12))
abline(v = med_neff, col = BLUE, lwd = 1.6)
text(med_neff + 0.15, max(h$counts) * 1.04,
     sprintf("median %.2f of 8", med_neff), adj = 0, col = BLUE, cex = 0.82)
mtext("A", side = 3, line = 0.4, adj = -0.16, font = 2, cex = 1.05)

## ---- B: MDE -------------------------------------------------------------
pan()
lx <- log10(ok$mde_fold)
h <- hist(lx, breaks = 60, plot = FALSE)
plot(h, col = "#e8d3d1", border = NA, main = "", xlab = "Minimum detectable fold change",
     ylab = "Genes", axes = FALSE, xlim = log10(c(1, 20)),
     ylim = c(0, max(h$counts) * 1.16))
at <- c(1, 1.5, 2, 3, 5, 10, 20)
axis(1, at = log10(at), labels = at)
axis(2)
abline(v = log10(med_mde), col = RED, lwd = 1.6)
abline(v = log10(mde_fdr), col = RED, lwd = 1.4, lty = 2)
text(log10(med_mde) - 0.03, max(h$counts) * 1.09, sprintf("%.2f", med_mde),
     adj = 1, col = RED, cex = 0.78)
text(log10(mde_fdr) + 0.03, max(h$counts) * 1.09, sprintf("%.2f (FDR)", mde_fdr),
     adj = 0, col = RED, cex = 0.78)
mtext("B", side = 3, line = 0.4, adj = -0.16, font = 2, cex = 1.05)

## ---- C: MDE vs depth ----------------------------------------------------
pan()
x <- dc$cells_per_donor_per_condition
plot(x, dc$median_mde_fold_all_input, log = "xy", type = "b", pch = 19, cex = 0.7,
     col = BLUE, xlab = "Cells per donor and condition",
     ylab = "Median MDE (fold change)", xaxt = "n",
     ylim = c(1.4, 14), xlim = c(8, 3200))
axis(1, at = c(10, 25, 50, 100, 250, 500, 1000, 2500), labels = c(10, 25, 50, 100, 250, 500, 1000, 2500))
points(x[length(x)], dc$median_mde_fold_all_input[length(x)], pch = 21, bg = RED,
       col = RED, cex = 1.5, lwd = 1.2)
text(2300, 7.2, sprintf("observed depth\n%.2f-fold", dc$median_mde_fold_all_input[length(x)]),
     adj = 1, col = RED, cex = 0.76)
mtext("C", side = 3, line = 0.4, adj = -0.16, font = 2, cex = 1.05)

## ---- D: n_eff vs depth --------------------------------------------------
pan()
plot(x, dc$median_n_eff_ok, log = "x", type = "b", pch = 19, cex = 0.7, col = GREY,
     xlab = "Cells per donor and condition", ylab = "Median effective donors",
     xaxt = "n", ylim = c(3.4, 6.2), xlim = c(8, 3200))
axis(1, at = c(10, 25, 50, 100, 250, 500, 1000, 2500), labels = c(10, 25, 50, 100, 250, 500, 1000, 2500))
points(x[length(x)], dc$median_n_eff_ok[length(x)], pch = 21, bg = RED, col = RED,
       cex = 1.5, lwd = 1.2)
mtext("D", side = 3, line = 0.4, adj = -0.16, font = 2, cex = 1.05)

dev.off()
cat("figure written:", file.path(OUT, "Figure1_v200.tif"), "\n")
cat("median n_eff", med_neff, "median MDE", med_mde, "\n")
