#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 08_hierarchicell_compare.R  (cloud host)
# Same dataset, same quantity: power to detect a two-fold change with 8 donors
# per group and 300 cells per donor. Hierarchicell simulates from parameters it
# estimates from GSE96583 and tests with a MAST hurdle model; effn's number for
# the same design is in out/unpaired_check.txt.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(hierarchicell)
})

OUT <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}
set.seed(20261021)

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
say("loaded ", nrow(M), " genes x ", ncol(M), " cells")

## top genes by total counts, then a fixed number of cells per donor
tot <- Matrix::rowSums(M)
genes <- names(sort(tot, decreasing = TRUE))[seq_len(300)]
keep <- unlist(lapply(split(seq_along(ob$donor), ob$donor), function(ix) {
  ix[seq_len(min(300L, length(ix)))]
}))
sub <- M[genes, keep, drop = FALSE]
say("subset: ", nrow(sub), " genes x ", ncol(sub), " cells, ",
    length(unique(ob$donor[keep])), " donors")

df <- data.frame(cell = colnames(sub), sample = ob$donor[keep],
                 as.matrix(t(sub)), check.names = FALSE)
rm(sub)
say("data frame: ", nrow(df), " x ", ncol(df), "; mem ",
    round(sum(gc()[, 2]), 0), " MB")

expr <- filter_counts(df, gene_thresh = 0, cell_thresh = 0)
say("filtered: ", nrow(expr), " cells, ", ncol(expr) - 2L, " genes")
ds <- compute_data_summaries(expr, type = "Raw")
say("summaries done: ", paste(names(ds), collapse = ", "))
saveRDS(ds, file.path(OUT, "hc_data_summaries.rds"))

msgs <- character(0)
grab <- function(m) {
  txt <- conditionMessage(m)
  msgs <<- c(msgs, txt)
  invokeRestart("muffleMessage")
}

t0 <- Sys.time()
res <- withCallingHandlers(
  power_hierarchicell(ds, n_genes = 100, n_per_group = 8,
                      cells_per_case = 300, cells_per_control = 300,
                      pval = 0.05, foldchange = 2),
  message = grab
)
say("power_hierarchicell took ",
    round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")
saveRDS(res, file.path(OUT, "hc_power.rds"))
writeLines(msgs, file.path(OUT, "hc_power_messages.txt"))
cat("--- messages ---\n"); cat(msgs, sep = "\n"); cat("\n--- structure ---\n")
str(res, max.level = 1)
writeLines(capture.output(str(res, max.level = 2)),
           file.path(OUT, "hc_power_str.txt"))
say("HCDONE")
