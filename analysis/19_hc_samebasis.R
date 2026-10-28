#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 19_hc_samebasis.R  (cloud host)
#
# v1.1.2 里有一条硬伤：同一段文字里摆出 Hierarchicell 的 0.11 与 effn 的 0.99，
# 结论相反却不解释。两者问的根本不是同一个问题：
#   Hierarchicell : 给定一群基因、给定每组供体数与细胞数，MAST hurdle 检验的
#                   平均 power（它自己模拟基因、自己造效应）
#   effn          : 给定一个基因，在 80% power 下能看见的最小倍数
# 把两个数并排放，等于把"平均 power"和"单个基因的检测下限"当同一件事。
#
# 本脚本改成可比的形式：两边都回答"power 随供体数怎么变"，倍数固定 2 倍。
#   effn 侧：由 sigma_obs 与非中心 t 直接算 power(delta = log2(2), n)
#   HC   侧：power_hierarchicell(foldchange = 2, n_per_group = n)
# 若两条曲线同向上升，说明二者一致；水平差异归因于检验与基因模型不同，
# 而不是谁对谁错。这是能写进论文的表述。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix); library(effn); library(hierarchicell)
})
OUT <- Sys.getenv("EFFN_OUT", "/root/effn/out")
say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush(stdout()) }

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")

## ---- effn 侧：非中心 t 的 power ----------------------------------------
pg <- utils::read.delim(gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                        stringsAsFactors = FALSE)
tot <- Matrix::rowSums(M)
top300 <- names(sort(tot, decreasing = TRUE))[seq_len(300)]
sub <- pg[match(top300, pg$gene), ]
sub <- sub[!is.na(sub$sigma_obs) & sub$mu_cpm >= 0.1, ]
say("effn side: ", nrow(sub), " genes with usable sigma_obs")

power_at <- function(delta_log2, n, sigma) {
  ncp <- delta_log2 * sqrt(n) / sigma
  crit <- qt(0.975, n - 1L)
  pt(-crit, n - 1L, ncp = ncp) + pt(crit, n - 1L, ncp = ncp, lower.tail = FALSE)
}
nl <- c(8L, 16L, 32L)
ef <- data.frame(
  n_per_group = nl,
  power_2fold = vapply(nl, function(n)
    mean(power_at(log2(2), n, sub$sigma_obs), na.rm = TRUE), numeric(1))
)
print(ef, digits = 4)

## ---- Hierarchicell 侧 ---------------------------------------------------
ds_file <- file.path(OUT, "hc_data_summaries.rds")
if (!file.exists(ds_file)) {
  keep <- unlist(lapply(split(seq_along(ob$donor), ob$donor), function(ix)
    ix[seq_len(min(300L, length(ix)))]))
  genes <- names(sort(tot, decreasing = TRUE))[seq_len(300)]
  sm <- M[genes, keep, drop = FALSE]
  df <- data.frame(cell = colnames(sm), sample = ob$donor[keep],
                   as.matrix(t(sm)), check.names = FALSE)
  rm(sm); gc()
  expr <- filter_counts(df, gene_thresh = 0, cell_thresh = 0)
  ds <- compute_data_summaries(expr, type = "Raw")
  saveRDS(ds, ds_file)
  rm(df, expr); gc()
} else {
  ds <- readRDS(ds_file)
}
say("Hierarchicell summaries ready")

grab <- function(msgs) function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
pick <- function(msgs, lab) {
  v <- grep(lab, msgs, value = TRUE)
  if (!length(v)) return(NA_real_)
  as.numeric(sub("^.*: ", "", v[length(v)]))
}

rows <- list()
for (n in nl) {
  msgs <- character(0)
  ## 注意：handler 必须定义在循环作用域里，<<- 才能写回 msgs；
  ## 若写成 grab(msgs) 返回闭包，<<- 会落在闭包自己的调用帧里，msgs 永远是空的。
  h <- function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  t0 <- Sys.time()
  set.seed(20261021)
  invisible(withCallingHandlers(
    power_hierarchicell(ds, n_genes = 100, n_per_group = n,
                        cells_per_case = 300, cells_per_control = 300,
                        pval = 0.05, foldchange = 2),
    message = h))
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  rows[[length(rows) + 1L]] <- data.frame(
    n_per_group = n,
    hc_hurdle = pick(msgs, "Hurdle Power"),
    hc_continuous = pick(msgs, "Continuous Power"),
    hc_discrete = pick(msgs, "Discrete Power"),
    minutes = el)
  say("HC n=", n, " hurdle=", rows[[length(rows)]]$hc_hurdle, " (", el, " min)")
}
hc <- do.call(rbind, rows)
cmp <- merge(ef, hc, by = "n_per_group")
print(cmp, digits = 4)
utils::write.table(cmp, file.path(OUT, "t_hc_samebasis.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("HCSAMEDONE\n")
