## 22b_jackknife_export.R -- 导出包内 jackknife 的 8 个 delete-one 值
##
## 背景：effn 的 jack 分支（effn.R:362-390）不重跑整个估计，而是固定全量
##       sigma_tech（t_donor），只按 delete-one 后的 donor 级方差重算
##       sigma_bio、rho 与 MDE。正文 3.91（2.39-5.43）与 1.63（1.48-1.81）
##       正是这一实现的输出，Table S10 必须与它同源，否则表与正文又对不上。
##
## 输出：out/jack_export.json / out/jack_export.tsv

LIB <- "C:/rtmp/effnlib"
.libPaths(c(LIB, .libPaths()))
suppressPackageStartupMessages({ library(Matrix); library(effn) })

ROOT <- "C:/rtmp/effnwork"
OUT  <- file.path(ROOT, "out")
say  <- function(...) { cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = ""); flush(stdout()) }

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M  <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
donor <- ob$donor; grp <- ob$group
say("loaded ", nrow(M), " x ", ncol(M))

say("fitting with jack = TRUE (package implementation) ...")
t0 <- Sys.time()
r <- effn(M, group = grp, donor = donor, n_boot = 200L, n_tech = 40L,
          seed = 20261021L, jack = TRUE)
say("fit took ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")

s <- r$summary
say(sprintf("full: n_ok %d | rho %.6f | n_eff %.6f | MDE %.6f (%.6f-fold) | share %.6f",
            s$n_ok, s$median_rho_tech, s$median_n_eff, s$median_mde_log2,
            s$median_mde_fold, s$loss_depth_share))
say("stored CI n_eff : ", paste(sprintf("%.6f", s$median_n_eff_ci), collapse = " to "))
say("stored CI MDE   : ", paste(sprintf("%.6f", s$median_mde_log2_ci), collapse = " to "),
    "  (", paste(sprintf("%.4f", 2^s$median_mde_log2_ci), collapse = " to "), "-fold )")
say("ci_method       : ", s$ci_method)

jk <- r$jack
stopifnot(!is.null(jk))
say("delete-one values (n_eff): ", paste(sprintf("%.4f", jk$n_eff), collapse = " "))
say("delete-one values (mde)  : ", paste(sprintf("%.4f", jk$mde), collapse = " "))

donors <- sort(unique(donor))
tdf <- data.frame(donor = donors, n_eff = jk$n_eff, mde_log2 = jk$mde,
                  mde_fold = 2^jk$mde)
print(tdf, digits = 6)
write.table(tdf, file.path(OUT, "jack_export.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

## 与正文存储的区间对照：应当逐位相同
chk <- c(
  n_eff_lo = s$median_n_eff_ci[1], n_eff_hi = s$median_n_eff_ci[2],
  mde_lo   = s$median_mde_log2_ci[1], mde_hi = s$median_mde_log2_ci[2]
)
say("--- interval check ---")
say(sprintf("n_eff CI %.6f to %.6f  (stated in text: 2.39 to 5.43)", chk[1], chk[2]))
say(sprintf("MDE   CI %.6f to %.6f  (%.4f to %.4f fold; stated: 1.48 to 1.81)",
            chk[3], chk[4], 2^chk[3], 2^chk[4]))

lines <- c(
  "{",
  sprintf('  "n_donor": %d,', length(donors)),
  sprintf('  "full": {"n_ok": %d, "rho": %.10f, "n_eff": %.10f, "mde_log2": %.10f, "mde_fold": %.10f, "share": %.10f},',
          s$n_ok, s$median_rho_tech, s$median_n_eff, s$median_mde_log2, s$median_mde_fold, s$loss_depth_share),
  sprintf('  "ci": {"n_eff_lo": %.10f, "n_eff_hi": %.10f, "mde_log2_lo": %.10f, "mde_log2_hi": %.10f, "mde_fold_lo": %.10f, "mde_fold_hi": %.10f},',
          chk[1], chk[2], chk[3], chk[4], 2^chk[3], 2^chk[4]),
  sprintf('  "delete_one": {"n_eff": [%s], "mde_log2": [%s]},',
          paste(sprintf("%.10f", jk$n_eff), collapse = ", "),
          paste(sprintf("%.10f", jk$mde), collapse = ", ")),
  sprintf('  "spread": {"n_eff_min": %.10f, "n_eff_max": %.10f, "mde_log2_min": %.10f, "mde_log2_max": %.10f},',
          min(jk$n_eff), max(jk$n_eff), min(jk$mde), max(jk$mde)),
  sprintf('  "ci_method": "%s"', s$ci_method),
  "}"
)
writeLines(lines, file.path(OUT, "jack_export.json"))
say("JACK_EXPORT_OK")
