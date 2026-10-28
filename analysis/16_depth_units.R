#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 16_depth_units.R  (cloud host)
#
# 目的：把深度曲线的横轴单位改对，并把 MDE 一起画出来。
#
# v1.1.2 的问题（两处，互相牵连）：
#   1. 横轴 cc 来自 table(donor)，即"该供体两个条件的细胞总数"；但校准臂和
#      800 细胞子样本用的是"每条件细胞数"。同一个词 cells per donor 两个分母。
#      实测：两条件合计最小 1,042（供体 107），每条件最小 461（供体 1039 ctrl）。
#      所以深度图上的 "100 cells" 实际等于子样本分析里的 "50 cells per condition"。
#   2. 只输出 rho 的占比（depth_share），没有输出 MDE，导致"深度到底买不买得到
#      分辨力"这件事只能靠 rho 间接推断——而 rho 与 MDE 不是单调对应的。
#
# 本脚本在下采样时按 (donor, group) 分层，于是横轴的字面含义就是"每条件每供体
# 细胞数"，与校准臂、子样本分析一致。同时输出 MDE 与 n_eff。
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

OUT <- Sys.getenv("EFFN_OUT", "/root/effn/out")
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
donor <- ob$donor; grp <- ob$group
key <- paste(donor, grp, sep = "\r")
say("loaded ", nrow(M), " genes x ", ncol(M), " cells")

## 输入基因集：主分析里 status == ok 的那 6,458 个。注意这里只用于"报告哪些
## 基因"，不用于限制 effn 的输入矩阵 —— CPM 的分母是库大小，若先把矩阵裁到
## 6,458 个基因再算，库大小会变，同一批基因的 rho 与 MDE 都会平移（实测满深度
## 那一点会给出 1.82 而不是 1.63）。所以每个深度都喂全基因矩阵，只在汇总时
## 取 top 这批基因。
main <- utils::read.delim(gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                          stringsAsFactors = FALSE)
top <- main$gene[main$status == "ok"]
say("report set: ", length(top), " genes estimable at the observed depth")

cnt_pc <- table(key)
min_pc <- min(cnt_pc)
say("cells per (donor, condition): min ", min_pc, ", median ", median(cnt_pc),
    ", max ", max(cnt_pc))
say("cells per donor (both conditions): min ", min(table(donor)))

grid <- c(10, 25, 50, 100, 250, 500, 1000, 2500)
grid <- sort(unique(pmin(grid, min_pc)))
## 再加一个"不降采样"的锚点：曲线必须能落在真实深度上，否则最小供体（461）会
## 被当成满深度，读者无法把曲线和正文的 1.63 倍对上。
grid <- c(grid, NA_real_)
say("depth grid (cells per donor per condition): ",
    paste(grid, collapse = ", "), " (NA = all cells)")

rows <- list()
for (cc in grid) {
  set.seed(20261021L)
  keep <- if (is.na(cc)) seq_along(donor) else
    unlist(lapply(split(seq_along(donor), key), function(ii)
      sample(ii, min(cc, length(ii)))), use.names = FALSE)
  Ms <- M[, keep, drop = FALSE]                    # 全基因，保证 CPM 分母一致
  t0 <- Sys.time()
  r <- effn(Ms, group = grp[keep], donor = donor[keep],
            n_boot = 0L, n_tech = 40L, seed = 20261021L)
  s <- r$summary
  pg <- r$per_gene[match(top, r$per_gene$gene), ]   # 只在汇总时取报告集
  ## 同一批基因，三个口径并排：可估子集 / 全部报告集
  ok <- pg$status == "ok"
  ts <- qt(0.975, 7) + qt(0.80, 7)
  mde_obs <- ts * pg$sigma_obs / sqrt(8)          # 与 s2_bio 口径代数等价
  rows[[length(rows) + 1L]] <- data.frame(
    cells_per_donor_per_condition = cc,
    cells_per_donor_total = if (is.na(cc)) NA_real_ else 2 * cc,
    n_cells_used = length(keep),
    n_ok = sum(ok),
    median_rho_ok = stats::median(pg$rho_tech[ok], na.rm = TRUE),
    median_n_eff_ok = stats::median(pg$n_eff[ok], na.rm = TRUE),
    median_mde_ok = stats::median(pg$mde_log2[ok], na.rm = TRUE),
    median_mde_fold_ok = 2^stats::median(pg$mde_log2[ok], na.rm = TRUE),
    median_mde_all_input = stats::median(mde_obs, na.rm = TRUE),
    median_mde_fold_all_input = 2^stats::median(mde_obs, na.rm = TRUE),
    tech_share_ok = stats::median((pg$rho_tech / (1 + pg$rho_tech))[ok], na.rm = TRUE),
    seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  )
  d <- rows[[length(rows)]]
  say("depth ", cc, "/cond -> n_ok ", d$n_ok,
      " mde_ok ", round(d$median_mde_fold_ok, 3),
      " mde_all ", round(d$median_mde_fold_all_input, 3),
      " n_eff ", round(d$median_n_eff_ok, 2),
      " (", d$seconds, " s)")
}

dc <- do.call(rbind, rows)
utils::write.table(dc, file.path(OUT, "depth_curve_per_condition.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
print(dc, digits = 4)
cat("DEPTHDONE\n")
