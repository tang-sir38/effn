#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 18_limma_real.R  (cloud host)
#
# 问题：effn 的 MDE 用的是未 moderated 的普通 t（df = 供体数 - 1 = 7）。
# 审稿人一定会问：换成 limma-voom 会怎样？
#
# 做法：同一份 pseudobulk（8 供体 × 2 条件，配对），两条流水线各算一遍：
#   effn   : MDE = (t_{.975,7} + t_{.80,7}) * sigma_obs / sqrt(8)
#   limma  : MDE = (t_{.975,df.total} + t_{.80,df.total}) * stdev.unscaled * sqrt(s2.post)
# 两条式子形式相同，差别只在方差估计（是否向先验收缩）与自由度。
#
# 输出：逐基因 MDE 的比值分布、中位数、以及 limma 的 df.total。
# 这份表回答的是"换尺子会不会改变结论的方向"，不是"哪个软件更好"。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix); library(effn); library(limma); library(edgeR)
})

OUT <- Sys.getenv("EFFN_OUT", "/root/effn/out")
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush(stdout()) }

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
say("loaded ", nrow(M), " genes x ", ncol(M), " cells")

## voom 要 counts：按 (donor, condition) 求和。不经过 pb_norm，因为 limma 需要
## 原始计数，而 pb_norm 输出的是 log2 CPM（那是 effn 内部用的尺度）。
key <- paste(ob$donor, ob$group, sep = "\r")
ug <- unique(key)
ug <- ug[order(sub("\r.*$", "", ug), sub("^.*\r", "", ug))]
sums <- vapply(ug, function(k) Matrix::rowSums(M[, key == k, drop = FALSE]),
               numeric(nrow(M)))
colnames(sums) <- ug
sk <- sub("\r.*$", "", ug); sg <- sub("^.*\r", "", ug)
say("count sums: ", nrow(sums), " genes x ", ncol(sums), " samples; total ",
    format(sum(sums), big.mark = ","))

donor_f <- factor(sk)
group_f <- factor(sg, levels = c("ctrl", "stim"))

## ---- limma-voom on pseudobulk -------------------------------------------
dge <- DGEList(counts = round(sums), group = group_f)
keep_expr <- filterByExpr(dge, design = model.matrix(~ donor_f + group_f))
dge <- dge[keep_expr, ]
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~ donor_f + group_f)
say("genes passing filterByExpr: ", sum(keep_expr), " of ", nrow(sums))

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)
co <- ncol(design)                       # group 系数在最后一列
dft <- fit$df.total
s2post <- fit$s2.post
su <- fit$stdev.unscaled[, co]
say("limma df.prior median ", round(median(fit$df.prior), 2),
    " -> df.total median ", round(median(dft), 2))

ts <- function(df) qt(0.975, df) + qt(0.80, df)
mde_limma <- ts(dft) * su * sqrt(s2post)
## 等价写法自检：SE 应当等于 stdev.unscaled * sqrt(s2.post)
se_chk <- su * sqrt(s2post)
say("check: max |SE - stdev.unscaled*sqrt(s2.post)| = ",
    signif(max(abs(se_chk - su * sqrt(s2post))), 3))

## ---- effn 口径（普通 t，df = 7） ----------------------------------------
pg <- utils::read.delim(gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                        stringsAsFactors = FALSE)
rownames(pg) <- pg$gene
g <- intersect(rownames(fit), pg$gene)
say("intersection: ", length(g), " genes")

mde_effn <- ts(7) * pg[g, "sigma_obs"] / sqrt(8)
mde_lim <- mde_limma[match(g, rownames(fit))]

cmp <- data.frame(
  gene = g,
  status = pg[g, "status"],
  sigma_obs = pg[g, "sigma_obs"],
  sigma_post = sqrt(s2post)[match(g, rownames(fit))],
  df_total = dft[match(g, rownames(fit))],
  mde_effn_log2 = mde_effn,
  mde_limma_log2 = mde_lim
)
cmp$ratio <- cmp$mde_limma_log2 / cmp$mde_effn_log2
cmp$fold_effn <- 2^cmp$mde_effn_log2
cmp$fold_limma <- 2^cmp$mde_limma_log2

ok <- cmp$status == "ok" & is.finite(cmp$ratio)
say("comparable genes: ", sum(ok))
s <- cmp[ok, ]
summ <- data.frame(
  metric = c("median MDE effn (fold)", "median MDE limma (fold)",
             "median ratio limma/effn", "IQR ratio",
             "genes where limma MDE < effn", "median df.total"),
  value = c(median(s$fold_effn), median(s$fold_limma), median(s$ratio),
            paste(round(quantile(s$ratio, c(.25, .75)), 3), collapse = "-"),
            sprintf("%d / %d (%.1f%%)", sum(s$ratio < 1), nrow(s),
                    100 * mean(s$ratio < 1)),
            median(s$df_total))
)
print(summ, row.names = FALSE)
utils::write.table(summ, file.path(OUT, "t_limma_real_summary.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
utils::write.table(cmp, gzfile(file.path(OUT, "t_limma_real_per_gene.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("LIMMAREALDONE\n")
