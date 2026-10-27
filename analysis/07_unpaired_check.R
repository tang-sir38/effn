#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 07_unpaired_check.R  (cloud host)
# The same data, analysed as an unpaired case-control design. Purpose: price
# the pairing, and give Hierarchicell a like-for-like comparison (its simulator
# is case-control, not paired).
#
# The two runs do not keep the same genes: a gene can be estimable under one
# design and not the other, so comparing two summary medians compares two
# different gene sets as well as two designs. Every ratio below is therefore
# computed on the intersection of the two ok sets, and the size of that
# intersection is reported next to it.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

OUT <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}

ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
say("loaded ", nrow(M), " x ", ncol(M))

## Deliberate mis-specification, run to price it: treat each donor-group as an
## independent patient, which is what an analysis that ignores the pairing does.
## Unique labels keep the package's repeated-donor guard satisfied without
## changing any count.
donor_naive <- paste0(ob$donor, "_", ob$group)

t0 <- Sys.time()
rep_u <- effn(M, group = ob$group, donor = donor_naive, n_boot = 200L, n_tech = 40L,
              seed = 20261021L, paired = FALSE, verbose = TRUE)
say("unpaired effn took ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), " min")

saveRDS(rep_u, file.path(OUT, "gse96583_unpaired.rds"))
writeLines(jsonlite::toJSON(rep_u$summary, auto_unbox = TRUE, digits = 10, pretty = TRUE),
           file.path(OUT, "gse96583_unpaired_summary.json"))
utils::write.table(rep_u$per_gene, gzfile(file.path(OUT, "gse96583_unpaired_per_gene.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)

## ---- like-for-like: the same genes under both designs ---------------------
rep_p <- readRDS(file.path(OUT, "gse96583_report.rds"))
pp <- rep_p$per_gene; pu <- rep_u$per_gene
common <- intersect(pp$gene[pp$status == "ok"], pu$gene[pu$status == "ok"])
say("ok under both designs: ", length(common), " genes")

pu_c <- pu[match(common, pu$gene), ]
pp_c <- pp[match(common, pp$gene), ]
mde_p <- stats::median(pp_c$mde_log2, na.rm = TRUE)
mde_u <- stats::median(pu_c$mde_log2, na.rm = TRUE)
ratio_common <- mde_u / mde_p
say("median MDE on the common set: paired ", round(mde_p, 4),
    " unpaired ", round(mde_u, 4), " ratio ", round(ratio_common, 4))

## power at a two-fold change for the median gene, unpaired, n = 8 per group.
## The quantile uses the design df (2n - 2); the effective count sets the SE.
sb <- rep_u$summary$median_sigma_bio
rho <- rep_u$summary$median_rho_tech
ne <- 8 / (1 + rho)
df <- 2 * 8 - 2
pow <- function(delta, sb, ne, df, alpha = 0.05) {
  tc <- stats::qt(1 - alpha / 2, df)
  ncp <- delta / (sb * sqrt(2 / ne))
  stats::pt(-tc, df, ncp) + stats::pt(tc, df, ncp, lower.tail = FALSE)
}
p2 <- pow(log2(2), sb, ne, df)
p15 <- pow(log2(1.5), sb, ne, df)
cat(sprintf("median gene: sigma_bio=%.4f rho=%.4f n_eff=%.3f df=%d\n",
            sb, rho, ne, df))
cat(sprintf("power at 2.0-fold, n=8/group, unpaired: %.4f\n", p2))
cat(sprintf("power at 1.5-fold, n=8/group, unpaired: %.4f\n", p15))

writeLines(c(
  sprintf("n_common=%d", length(common)),
  sprintf("median_sigma_bio=%.6f", sb),
  sprintf("median_rho=%.6f", rho),
  sprintf("n_eff_at_8=%.6f", ne),
  sprintf("df_design=%d", df),
  sprintf("power_2fold_n8=%.6f", p2),
  sprintf("power_1.5fold_n8=%.6f", p15),
  sprintf("mde_unpaired_common=%.6f", mde_u),
  sprintf("mde_paired_common=%.6f", mde_p),
  sprintf("mde_ratio_unpaired_over_paired=%.6f", ratio_common),
  sprintf("mde_paired_allgenes=%.6f", rep_p$summary$median_mde_log2),
  sprintf("mde_unpaired_allgenes=%.6f", rep_u$summary$median_mde_log2)
), file.path(OUT, "unpaired_check.txt"))
say("UNPAIREDDONE")
