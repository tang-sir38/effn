#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 11_small_reproduce.R  (any machine, no network, no GSE96583 download)
#
# Reproduces the per-gene table of the 800-cell subsample that ships inside the
# package (inst/extdata/gse96583_small.rds) and writes the facts file the number
# ledger reads (03_RESULTS/small_dataset.txt).
#
# The md5 quoted in the paper is the md5 of the *distributed* file
# (small_per_gene.tsv.gz, shipped in the Supplementary Data). Re-running this
# script reproduces the medians exactly; the last printed digit of mde_fold can
# move on a different platform, because write.table prints doubles at 15
# significant digits and the two builds round a handful of them differently.
#
# Usage:
#   Rscript 11_small_reproduce.R [pkgdir] [outdir] [shipped_gz]
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args) >= 1) args[1] else "."
outdir <- if (length(args) >= 2) args[2] else tempdir()
shipped <- if (length(args) >= 3) args[3] else NA_character_

suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

rds <- file.path(pkg, "inst", "extdata", "gse96583_small.rds")
if (!file.exists(rds)) rds <- file.path(pkg, "effn", "inst", "extdata", "gse96583_small.rds")
if (!file.exists(rds)) rds <- system.file("extdata", "gse96583_small.rds", package = "effn")
if (!nzchar(rds) || !file.exists(rds)) {
  stop("cannot find gse96583_small.rds; pass the package directory as the first argument")
}
small <- readRDS(rds)
stopifnot(all(c("counts", "donor", "group") %in% names(small)))

tb <- table(small$donor, small$group)
cat("object:", nrow(small$counts), "genes x", ncol(small$counts), "cells;",
    length(unique(small$donor)), "donors\n")
stopifnot(ncol(tb) == 2L, all(tb == 50L))

t0 <- Sys.time()
rep_s <- effn(small$counts, group = small$group, donor = small$donor,
              n_boot = 200L, n_tech = 40L, seed = 20261021L)
el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

tsv <- file.path(outdir, "small_per_gene.tsv")
con <- file(tsv, "wb")            # binary: no CRLF translation, same bytes anywhere
write.table(rep_s$per_gene, con, sep = "\t", quote = FALSE, row.names = FALSE)
close(con)
md5_tsv <- unname(tools::md5sum(tsv))

cat(sprintf("elapsed: %.1f s\n", el))
cat("md5(per_gene.tsv, plain) =", md5_tsv, "\n")
cat(sprintf("median n_eff = %.10f\n", rep_s$summary$median_n_eff))
cat(sprintf("median mde_log2 = %.10f\n", rep_s$summary$median_mde_log2))
cat(sprintf("median rho_tech = %.10f\n", rep_s$summary$median_rho_tech))
cat(sprintf("n_ok = %d\n", rep_s$summary$n_ok))

## facts file, one key=value per line, read by 02_SCRIPTS/build_ledger.py
facts <- c(
  genes = nrow(small$counts),
  cells = ncol(small$counts),
  donors = length(unique(small$donor)),
  cells_per_donor_group = as.integer(tb[1, 1]),
  n_ok = rep_s$summary$n_ok,
  median_n_eff = sprintf("%.10f", rep_s$summary$median_n_eff),
  median_mde_log2 = sprintf("%.10f", rep_s$summary$median_mde_log2),
  median_rho_tech = sprintf("%.10f", rep_s$summary$median_rho_tech),
  md5_per_gene_tsv_plain = md5_tsv,
  secs = sprintf("%.1f", el),
  r_version = R.version.string,
  md5_per_gene_gz = if (is.na(shipped)) "NA" else unname(tools::md5sum(shipped))
)
writeLines(paste(names(facts), facts, sep = "="),
           file.path(outdir, "small_dataset.txt"))
cat("wrote", file.path(outdir, "small_dataset.txt"), "\n")
cat("SMALLREPRODONE\n")
