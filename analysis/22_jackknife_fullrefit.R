## 22_jackknife_gse96583.R -- 真 delete-one-donor，重建 Table S10
##
## 动机：v1.2.0 的 Table S10 标题写 "leave-one-donor-out fits"，但表内数字
##       （rho 0.601 / n_eff 4.997 / n_ok 72）来自 07_rho_uncertainty.R 的
##       30 次独立模拟，与正文 GSE96583 的 3.91 / 0.71 / 6,458 不是同一套。
##       本脚本对真实数据做 8 次拟合并给出 jackknife 区间。
##
## 口径：与 02_effn_main.R 一致（paired, n_tech = 40, seed = 20261021），
##       仅把 n_boot 设为 0（jackknife 不需要 bootstrap）。
##
## 输出：out/jackknife_full.tsv        全量拟合的摘要
##       out/jackknife_runs.tsv        8 次 delete-one 的逐次摘要
##       out/jackknife_summary.json    jackknife 点估计、SE、95% CI
##       out/jackknife_per_gene.tsv.gz 每次拟合的逐基因估计（用于核对）

LIB <- "C:/rtmp/effnlib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({ library(Matrix); library(effn) })

ROOT <- "C:/rtmp/effnwork"
DAT  <- file.path(ROOT, "data")
OUT  <- file.path(ROOT, "out")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

say <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush(stdout()) }

## ---- 1. 载入重建的 counts ------------------------------------------------
ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M  <- methods::as(ob$counts, "CsparseMatrix"); M <- as(M, "dgCMatrix")
donor <- ob$donor; grp <- ob$group
say("loaded ", nrow(M), " x ", ncol(M), " | donors ", length(unique(donor)))

## jsonlite 不可用，改用极简解析（文件是本包自己写的平铺 JSON）
jget <- function(txt, key) {
  m <- regmatches(txt, regexpr(paste0('"', key, '"\\s*:\\s*\\[?[-0-9.eE+]+'), txt))
  as.numeric(sub(paste0('.*:\\s*\\[?'), '', m))
}
jtxt <- paste(readLines(file.path(DAT, "gse96583_summary.json"), warn = FALSE), collapse = " ")
ref <- list(
  n_ok              = jget(jtxt, "n_ok"),
  median_rho_tech   = jget(jtxt, "median_rho_tech"),
  median_n_eff      = jget(jtxt, "median_n_eff"),
  median_mde_log2   = jget(jtxt, "median_mde_log2"),
  median_mde_fold   = jget(jtxt, "median_mde_fold"),
  median_n_eff_ci   = c(jget(jtxt, "median_n_eff") , NA),
  median_mde_log2_ci= c(NA, NA)
)
## 区间单独取（数组形式）
ci <- regmatches(jtxt, regexpr('"median_n_eff_ci"\\s*:\\s*\\[[^]]+\\]', jtxt))
ref$median_n_eff_ci <- if (length(ci)) as.numeric(strsplit(gsub('.*\\[|\\]', '', ci), ",")[[1]]) else c(NA, NA)
ci2 <- regmatches(jtxt, regexpr('"median_mde_log2_ci"\\s*:\\s*\\[[^]]+\\]', jtxt))
ref$median_mde_log2_ci <- if (length(ci2)) as.numeric(strsplit(gsub('.*\\[|\\]', '', ci2), ",")[[1]]) else c(NA, NA)
say("reference: n_eff ", ref$median_n_eff, " | MDE log2 ", ref$median_mde_log2,
    " | rho ", ref$median_rho_tech, " | n_ok ", ref$n_ok)

## effn() 需要矩阵与标签对齐；这里用列索引传入子集
run_effn <- function(cols, tag) {
  t0 <- Sys.time()
  r <- try(effn(M[, cols, drop = FALSE], group = grp[cols], donor = donor[cols],
                n_boot = 0L, n_tech = 40L, seed = 20261021L), silent = TRUE)
  if (inherits(r, "try-error")) {
    say(tag, " FAILED: ", as.character(r)); return(NULL)
  }
  s <- r$summary
  say(sprintf("%s done in %.1f min | n_ok %d | rho %.4f | n_eff %.4f | MDE %.4f (%.4f-fold)",
              tag, as.numeric(difftime(Sys.time(), t0, units = "mins")),
              s$n_ok, s$median_rho_tech, s$median_n_eff, s$median_mde_log2, s$median_mde_fold))
  list(summary = s, per_gene = r$per_gene)
}

## ---- 2. 全量口径核对 -----------------------------------------------------
say("--- full fit (all 8 donors) ---")
full <- run_effn(seq_len(ncol(M)), "full")
if (is.null(full)) stop("full fit failed")

chk <- data.frame(
  quantity = c("n_ok", "median_rho_tech", "median_n_eff", "median_mde_log2", "median_mde_fold"),
  rebuilt  = c(full$summary$n_ok, full$summary$median_rho_tech, full$summary$median_n_eff,
               full$summary$median_mde_log2, full$summary$median_mde_fold),
  recorded = c(ref$n_ok, ref$median_rho_tech, ref$median_n_eff,
               ref$median_mde_log2, ref$median_mde_fold)
)
chk$abs_diff <- abs(chk$rebuilt - chk$recorded)
chk$match    <- chk$abs_diff < 1e-6
print(chk, digits = 8)
write.table(chk, file.path(OUT, "jackknife_full_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!all(chk$match)) {
  say("WARNING: rebuilt counts do not reproduce the recorded summary to 1e-6.")
  say("         Jackknife will still run, but the discrepancy must be resolved before use.")
} else {
  say("full fit reproduces the recorded summary exactly")
}

## ---- 3. 8 次 delete-one-donor -------------------------------------------
donors <- sort(unique(donor))
say("--- jackknife over ", length(donors), " donors ---")

runs <- list(); pgs <- list()
for (k in donors) {
  idx <- which(donor != k)
  r <- run_effn(idx, paste0("drop ", k))
  if (is.null(r)) next
  runs[[as.character(k)]] <- data.frame(
    dropped = k, n_donor_left = length(unique(donor[idx])),
    n_cells = length(idx), n_ok = r$summary$n_ok,
    rho = r$summary$median_rho_tech, share = r$summary$loss_depth_share,
    n_eff = r$summary$median_n_eff, mde_log2 = r$summary$median_mde_log2,
    mde_fold = r$summary$median_mde_fold
  )
  pg <- r$per_gene
  pg$dropped <- k
  pgs[[as.character(k)]] <- pg[, c("dropped", "gene", "n_eff", "rho_tech", "mde_log2", "status")]
}

R <- do.call(rbind, runs)
rownames(R) <- NULL
print(R, digits = 5)
write.table(R, file.path(OUT, "jackknife_runs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

pgall <- do.call(rbind, pgs)
gz <- gzfile(file.path(OUT, "jackknife_per_gene.tsv.gz"), "w")
write.table(pgall, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)
say("per-gene estimates written for ", nrow(pgall), " rows")

## ---- 4. jackknife 点估计与区间 ------------------------------------------
JK <- function(x) {
  n <- length(x); m <- mean(x)
  se <- sqrt((n - 1) / n * sum((x - m)^2))
  c(mean = m, se = se, lo = m - qt(0.975, df = n - 1) * se, hi = m + qt(0.975, df = n - 1) * se)
}
jk_neff <- JK(R$n_eff)
jk_mde  <- JK(R$mde_log2)
jk_rho  <- JK(R$rho)
jk_nok  <- JK(R$n_ok)

jk_lines <- c(
  "{",
  sprintf('  "method": "jackknife (delete one donor)",'),
  sprintf('  "n_donors": %d,', length(donors)),
  sprintf('  "n_fits": %d,', nrow(R)),
  sprintf('  "full": {"n_ok": %d, "rho": %.10f, "n_eff": %.10f, "mde_log2": %.10f, "mde_fold": %.10f, "share": %.10f},',
          full$summary$n_ok, full$summary$median_rho_tech, full$summary$median_n_eff,
          full$summary$median_mde_log2, full$summary$median_mde_fold, full$summary$loss_depth_share),
  sprintf('  "jackknife": {'),
  sprintf('    "n_eff": {"mean": %.10f, "se": %.10f, "lo": %.10f, "hi": %.10f},', jk_neff[1], jk_neff[2], jk_neff[3], jk_neff[4]),
  sprintf('    "mde_log2": {"mean": %.10f, "se": %.10f, "lo": %.10f, "hi": %.10f},', jk_mde[1], jk_mde[2], jk_mde[3], jk_mde[4]),
  sprintf('    "mde_fold": {"mean": %.10f, "lo": %.10f, "hi": %.10f},', 2^jk_mde[1], 2^jk_mde[3], 2^jk_mde[4]),
  sprintf('    "rho": {"mean": %.10f, "se": %.10f, "lo": %.10f, "hi": %.10f},', jk_rho[1], jk_rho[2], jk_rho[3], jk_rho[4]),
  sprintf('    "n_ok": {"mean": %.4f, "se": %.4f, "lo": %.4f, "hi": %.4f}', jk_nok[1], jk_nok[2], jk_nok[3], jk_nok[4]),
  sprintf('  },'),
  sprintf('  "spread": {"n_eff_min": %.10f, "n_eff_max": %.10f, "mde_fold_min": %.10f, "mde_fold_max": %.10f, "rho_min": %.10f, "rho_max": %.10f},',
          min(R$n_eff), max(R$n_eff), min(R$mde_fold), max(R$mde_fold), min(R$rho), max(R$rho)),
  sprintf('  "ref_summary_ci": {"n_eff_lo": %.10f, "n_eff_hi": %.10f, "mde_log2_lo": %.10f, "mde_log2_hi": %.10f}', 
          ref$median_n_eff_ci[1], ref$median_n_eff_ci[2], ref$median_mde_log2_ci[1], ref$median_mde_log2_ci[2]),
  "}"
)
writeLines(jk_lines, file.path(OUT, "jackknife_summary.json"))

out <- list(
  method = "jackknife (delete one donor)",
  n_donors = length(donors),
  n_fits = nrow(R),
  full = list(n_ok = full$summary$n_ok, rho = full$summary$median_rho_tech,
              n_eff = full$summary$median_n_eff, mde_log2 = full$summary$median_mde_log2,
              mde_fold = full$summary$median_mde_fold,
              share = full$summary$loss_depth_share),
  jackknife = list(
    n_eff = as.list(jk_neff), mde_log2 = as.list(jk_mde),
    mde_fold = as.list(2^jk_mde), rho = as.list(jk_rho), n_ok = as.list(jk_nok)
  ),
  spread = list(n_eff_min = min(R$n_eff), n_eff_max = max(R$n_eff),
                mde_fold_min = min(R$mde_fold), mde_fold_max = max(R$mde_fold),
                rho_min = min(R$rho), rho_max = max(R$rho)),
  ref_summary_ci = list(n_eff_lo = ref$median_n_eff_ci[1], n_eff_hi = ref$median_n_eff_ci[2],
                        mde_log2_lo = ref$median_mde_log2_ci[1], mde_log2_hi = ref$median_mde_log2_ci[2])
)


say("--- jackknife summary ---")
say(sprintf("n_eff  : %.4f  SE %.4f  95%% CI %.4f to %.4f", jk_neff[1], jk_neff[2], jk_neff[3], jk_neff[4]))
say(sprintf("MDE    : %.4f  SE %.4f  95%% CI %.4f to %.4f  (%.4f to %.4f fold)",
            jk_mde[1], jk_mde[2], jk_mde[3], jk_mde[4], 2^jk_mde[3], 2^jk_mde[4]))
say(sprintf("rho    : %.4f  SE %.4f  95%% CI %.4f to %.4f", jk_rho[1], jk_rho[2], jk_rho[3], jk_rho[4]))
say(sprintf("n_ok   : %.2f  SE %.2f  95%% CI %.2f to %.2f", jk_nok[1], jk_nok[2], jk_nok[3], jk_nok[4]))
say("stored CI in v1.2.0 : n_eff ", ref$median_n_eff_ci, " | MDE log2 ", ref$median_mde_log2_ci)
say("JACKKNIFE_OK")
