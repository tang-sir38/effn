#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 02_effn_main.R  (cloud host)
# Stage B: run effn() on the assembled GSE96583 singlets, then re-run it with
# the same seed to prove byte-identical output. Resumable: if the summary file
# already exists the stage is skipped.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

OUT <- "/root/effn/out"
RDS <- file.path(OUT, "gse96583_singlets.rds")
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}

say("R ", as.character(getRversion()), " effn ", as.character(utils::packageVersion("effn")))
ob <- readRDS(RDS)
M <- methods::as(ob$counts, "CsparseMatrix")  # readMM gives dgTMatrix
M <- as(M, "dgCMatrix")
donor <- ob$donor; grp <- ob$group
say("loaded ", nrow(M), " genes x ", ncol(M), " cells; ", length(unique(donor)), " donors")
say("mem after load (MB): ", round(sum(gc()[, 2]), 0))

if (!file.exists(file.path(OUT, "gse96583_summary.json"))) {
  say("stage B1: effn ...")
  t0 <- Sys.time()
  rep1 <- effn(M, group = grp, donor = donor, n_boot = 200L, n_tech = 40L,
               seed = 20261021L, verbose = TRUE)
  say("effn took ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")
  saveRDS(rep1, file.path(OUT, "gse96583_report.rds"))
  utils::write.table(rep1$per_gene, gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(jsonlite::toJSON(rep1$summary, auto_unbox = TRUE, digits = 10, pretty = TRUE),
             file.path(OUT, "gse96583_summary.json"))
  capture.output(print(rep1), file = file.path(OUT, "gse96583_print.txt"))
  say("stage B1 done")
} else {
  say("stage B1 skipped: summary exists")
  rep1 <- readRDS(file.path(OUT, "gse96583_report.rds"))
}

if (!file.exists(file.path(OUT, "determinism.txt"))) {
  say("stage B2: determinism re-run ...")
  rep2 <- effn(M, group = grp, donor = donor, n_boot = 200L, n_tech = 40L,
               seed = 20261021L)
  same_pg <- identical(rep1$per_gene, rep2$per_gene)
  same_sm <- identical(rep1$summary, rep2$summary)
  same_bt <- identical(rep1$boot, rep2$boot)
  say("per_gene identical: ", same_pg, " summary identical: ", same_sm,
      " boot identical: ", same_bt)
  writeLines(c(paste0("per_gene_identical=", same_pg),
               paste0("summary_identical=", same_sm),
               paste0("boot_identical=", same_bt),
               paste0("identical=", same_pg && same_sm && same_bt)),
             file.path(OUT, "determinism.txt"))
  say("stage B2 done")
} else {
  say("stage B2 skipped: determinism.txt exists")
}

say("MAINDONE")
