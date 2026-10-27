#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 05_subsample_recheck.R
# Rebuilds the count matrix from the exported triplets and re-runs effn() with a
# fixed seed. Run on both machines; the printed hash of the per-gene table is
# what the manuscript compares.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

dir <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
tri <- utils::read.delim(gzfile(file.path(dir, "subsample_triplets.tsv.gz")),
                         header = FALSE)
genes <- readLines(gzfile(file.path(dir, "subsample_genes.txt.gz")))
cells <- utils::read.delim(gzfile(file.path(dir, "subsample_cells.tsv.gz")),
                           stringsAsFactors = FALSE)
M <- Matrix::sparseMatrix(i = tri[[1]], j = tri[[2]], x = tri[[3]],
                          dims = c(length(genes), nrow(cells)))
rownames(M) <- genes
colnames(M) <- cells$barcode
cat("matrix:", nrow(M), "x", ncol(M), " nnz", length(M@x), "\n")

t0 <- Sys.time()
r <- effn(M, group = cells$group, donor = cells$donor,
          n_boot = 50L, n_tech = 20L, seed = 20261021L)
cat("effn took", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")
print(r)

f <- file.path(dir, "subsample_per_gene.tsv")
utils::write.table(r$per_gene, f, sep = "\t", quote = FALSE, row.names = FALSE)
h <- tryCatch(tools::md5sum(f)[[1]], error = function(e) NA_character_)
cat("PER_GENE_MD5", h, "\n")
cat("MEDIAN_N_EFF", sprintf("%.10f", r$summary$median_n_eff),
    "MEDIAN_MDE", sprintf("%.10f", r$summary$median_mde_log2),
    "MEDIAN_RHO", sprintf("%.10f", r$summary$median_rho_tech), "\n")
cat("RECHECKDONE\n")
