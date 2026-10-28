#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 10_small_dataset.R  (cloud host)
# Builds the 800-cell subsample that ships with the package: 100 cells per
# donor, taken in column order so the subset is deterministic. The md5 of the
# per-gene table computed from this object is printed for the Supplementary.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})
OUT <- "/root/effn/out"
ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")

## 50 cells per donor AND condition, so the subset keeps both groups
key <- paste(ob$donor, ob$group, sep = "|")
idx <- unlist(lapply(split(seq_along(key), key),
                     function(ix) ix[seq_len(min(50L, length(ix)))]))
idx <- sort(idx)
small <- list(
  counts = M[, idx, drop = FALSE],
  donor = ob$donor[idx],
  group = ob$group[idx]
)
cat("subset:", nrow(small$counts), "genes x", ncol(small$counts), "cells;",
    length(unique(small$donor)), "donors\n")
tb <- table(small$donor, small$group)
stopifnot(ncol(tb) == 2L, all(tb == 50L))
cat("cells per donor x group: ", paste(unique(tb), collapse = "/"), "\n", sep = "")

saveRDS(small, file.path(OUT, "gse96583_small.rds"), compress = "xz")
cat("size KB:", round(file.size(file.path(OUT, "gse96583_small.rds")) / 1024, 1), "\n")

## the md5 the manuscript promises
t0 <- Sys.time()
rep_s <- effn(small$counts, group = small$group, donor = small$donor,
              n_boot = 200L, n_tech = 40L, seed = 20261021L)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
pg <- rep_s$per_gene
con <- gzfile(file.path(OUT, "small_per_gene.tsv.gz"), "w")
utils::write.table(pg, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
md5 <- tools::md5sum(file.path(OUT, "small_per_gene.tsv.gz"))
cat(sprintf("80d cells: %.1f s, md5(per_gene)=%s\n", el, unname(md5)))
cat(sprintf("median n_eff=%.10f mde=%.10f rho=%.10f\n",
            rep_s$summary$median_n_eff, rep_s$summary$median_mde_log2,
            rep_s$summary$median_rho_tech))
cat("SMALLDONE\n")
