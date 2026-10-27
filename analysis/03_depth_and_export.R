#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 03_depth_and_export.R  (cloud host)
# Stage C: depth curve by empirical downsampling (the title decision), plus the
# export of a subsample that the local machine can re-check independently.
#
# What the gene set is. Every depth is run on the SAME input set: the genes that
# are estimable at the observed depth. That set is fixed before any depth is
# looked at, so the input never depends on the depth.
#
# What the median is taken over. Estimability itself changes with depth, so a
# median taken over "whatever is still estimable at this depth" is a median over
# survivors, and survivors are the genes least sensitive to depth. That
# understates the penalty. Three quantities are therefore reported side by side:
#   core_median_*        median over survivors at that depth (the old number,
#                        kept because earlier versions quoted it)
#   fixed_median_*       median over the intersection: genes estimable at EVERY
#                        depth in the grid, so the same genes at every depth
#   core_n_ok/frac_lost  how much of the fixed input survives, so a reader can
#                        see the attrition that the survivor median hides
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
M <- methods::as(ob$counts, "CsparseMatrix")  # readMM gives dgTMatrix
M <- as(M, "dgCMatrix")
donor <- ob$donor; grp <- ob$group
say("loaded ", nrow(M), " genes x ", ncol(M), " cells")

## ---------------------------------------------------------------------------
if (!file.exists(file.path(OUT, "depth_curve.tsv"))) {
  say("stage C1: depth curve ...")
  cells_per <- table(donor)
  grid <- c(10, 25, 50, 100, 250, 500, 1000, 2500)
  grid <- sort(unique(pmin(grid, min(cells_per))))
  say("depth grid: ", paste(grid, collapse = ", "), " (smallest donor has ",
      min(cells_per), ")")

  ## fixed input set: the genes estimable at the observed depth, taken from the
  ## main report rather than from any property of a subsample
  main <- utils::read.delim(gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                            stringsAsFactors = FALSE)
  top <- main$gene[main$status == "ok"]
  say("fixed input set: ", length(top), " genes estimable at the observed depth")

  pg <- list(); meta <- list()
  for (cc in grid) {
    set.seed(20261021L)
    keep <- unlist(lapply(split(seq_along(donor), donor), function(ii)
      sample(ii, min(cc, length(ii)))), use.names = FALSE)
    Ms <- M[top, keep, drop = FALSE]
    t0 <- Sys.time()
    r <- effn(Ms, group = grp[keep], donor = donor[keep],
              n_boot = 0L, n_tech = 20L, seed = 20261021L)
    pg[[as.character(cc)]] <- r$per_gene[, c("gene", "status", "sigma_obs",
                                             "sigma_tech", "sigma_bio",
                                             "rho_tech", "n_eff", "mde_log2")]
    s <- r$summary
    meta[[as.character(cc)]] <- data.frame(
      cells_per_donor = cc, n_cells_used = length(keep),
      n_ok = s$n_ok,
      pct_tech_dominated = 100 * s$n_unestimable[["tech_dominated"]] / s$n_genes_tested,
      pct_low_resolution = 100 * s$n_unestimable[["low_resolution"]] / s$n_genes_tested,
      pct_low_expression = 100 * s$n_unestimable[["low_expression"]] / s$n_genes_tested,
      all_gene_median_rho = s$median_rho_tech,
      all_gene_median_n_eff = s$median_n_eff,
      all_gene_median_mde = s$median_mde_log2,
      seconds = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    )
    say("depth ", cc, " -> n_ok ", s$n_ok, " rho ", round(s$median_rho_tech, 3),
        " n_eff ", round(s$median_n_eff, 2), " mde ",
        round(s$median_mde_log2, 3), " (", meta[[as.character(cc)]]$seconds, " s)")
  }

  ## the intersection: genes marked ok at every depth in the grid
  okmat <- sapply(pg, function(d) d$status[match(top, d$gene)] == "ok")
  fixed <- top[rowSums(okmat) == ncol(okmat)]
  say("intersection set (ok at every depth): ", length(fixed), " genes")

  rows <- list()
  for (cc in grid) {
    d <- pg[[as.character(cc)]]
    d <- d[match(top, d$gene), ]
    okc <- d$status == "ok"
    fixed_ok <- okc & (top %in% fixed)
    nreq <- effn:::prescribe_n_vec(1, d$sigma_bio, d$rho_tech, paired = TRUE)
    rows[[length(rows) + 1L]] <- data.frame(
      cells_per_donor = cc,
      n_core = length(top),
      core_n_ok = sum(okc),
      core_median_sigma_tech = median(d$sigma_tech, na.rm = TRUE),
      core_median_sigma_bio = median(d$sigma_bio[okc], na.rm = TRUE),
      core_median_rho = median(d$rho_tech[okc], na.rm = TRUE),
      core_depth_share = median(d$rho_tech[okc] / (1 + d$rho_tech[okc]), na.rm = TRUE),
      core_median_n_eff = median(d$n_eff[okc], na.rm = TRUE),
      core_median_mde_log2 = median(d$mde_log2[okc], na.rm = TRUE),
      core_median_n_required_2x = median(nreq, na.rm = TRUE),
      core_frac_lost = 1 - mean(okc),
      fixed_n = sum(fixed_ok),
      fixed_median_rho = median(d$rho_tech[fixed_ok], na.rm = TRUE),
      fixed_depth_share = median(d$rho_tech[fixed_ok] / (1 + d$rho_tech[fixed_ok]),
                                 na.rm = TRUE),
      fixed_median_n_eff = median(d$n_eff[fixed_ok], na.rm = TRUE),
      fixed_median_mde_log2 = median(d$mde_log2[fixed_ok], na.rm = TRUE),
      fixed_median_n_required_2x = median(nreq[fixed_ok], na.rm = TRUE)
    )
  }
  dc <- cbind(do.call(rbind, rows), do.call(rbind, meta)[, -1, drop = FALSE])
  utils::write.table(dc, file.path(OUT, "depth_curve.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  print(dc, digits = 4)
  saveRDS(pg, file.path(OUT, "depth_per_gene.rds"))
  writeLines(top, file.path(OUT, "depth_core_genes.txt"))
  writeLines(fixed, file.path(OUT, "depth_fixed_genes.txt"))
  say("stage C1 done")
} else {
  say("stage C1 skipped")
}

## ---------------------------------------------------------------------------
## Exchange rate: what one donor is worth, in cells.
##
## Comparing medians across depths compares gene sets, because estimability
## changes with depth. The ratio below is therefore taken WITHIN gene: for every
## gene with a finite MDE at both depths, the factor by which its own MDE grows
## when the donor is thinned to d cells. The median of those per-gene factors is
## free of any composition effect. A design with n donors at depth d then has
## the same median MDE as a design with n * r^2 donors at full depth, where r is
## the median factor: the price of shallow sequencing expressed in donors.
if (!file.exists(file.path(OUT, "exchange_rate.tsv"))) {
  pg <- readRDS(file.path(OUT, "depth_per_gene.rds"))
  grid <- as.numeric(names(pg))
  ref <- pg[[as.character(max(grid))]]
  rows <- list()
  for (cc in sort(grid)) {
    d <- pg[[as.character(cc)]]
    m <- merge(d[, c("gene", "status", "mde_log2")],
               ref[, c("gene", "status", "mde_log2")], by = "gene",
               suffixes = c("", "_ref"))
    both <- m$status == "ok" & m$status_ref == "ok" &
      is.finite(m$mde_log2) & is.finite(m$mde_log2_ref)
    if (sum(both) < 10L) next
    r <- 2^(m$mde_log2[both] - m$mde_log2_ref[both])
    q <- stats::quantile(r, c(0.25, 0.5, 0.75), names = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(
      cells_per_donor = cc,
      n_genes_both = sum(both),
      median_fold_vs_full = q[2],
      iqr_low = q[1], iqr_high = q[3],
      ## a donor at this depth is worth 1 / r^2 of a full-depth donor, because
      ## the MDE falls as the square root of the donor count
      donor_equivalent_of_8 = 8 / q[2]^2,
      cells_per_full_depth_donor = cc * q[2]^2
    )
  }
  ex <- do.call(rbind, rows)
  utils::write.table(ex, file.path(OUT, "exchange_rate.tsv"), sep = "\t",
                     quote = FALSE, row.names = FALSE)
  print(ex, digits = 4)
  say("exchange rate written (within-gene ratios)")
} else {
  say("exchange rate skipped")
}

## ---------------------------------------------------------------------------
if (!file.exists(file.path(OUT, "subsample_cells.tsv.gz"))) {
  say("stage C2: subsample export ...")
  set.seed(20261022L)
  sub <- unlist(lapply(split(seq_along(donor), donor), function(ii)
    sample(ii, min(100, length(ii)))), use.names = FALSE)
  Ms <- M[, sub, drop = FALSE]
  tr <- summary(Ms)
  utils::write.table(data.frame(i = tr$i, j = tr$j, x = tr$x),
                     gzfile(file.path(OUT, "subsample_triplets.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeLines(rownames(Ms), gzfile(file.path(OUT, "subsample_genes.txt.gz")))
  utils::write.table(data.frame(barcode = ob$barcode[sub], donor = donor[sub],
                                group = grp[sub]),
                     gzfile(file.path(OUT, "subsample_cells.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  say("stage C2 done: ", nrow(Ms), " x ", ncol(Ms), " nnz ", length(Ms@x))
} else {
  say("stage C2 skipped")
}

if (!file.exists(file.path(OUT, "pseudobulk_log2cpm.tsv.gz"))) {
  say("stage C3: pseudobulk export ...")
  pb <- effn:::pb_norm(M, paste(donor, grp, sep = "|"))
  utils::write.table(as.data.frame(pb), gzfile(file.path(OUT, "pseudobulk_log2cpm.tsv.gz")),
                     sep = "\t", quote = FALSE, row.names = TRUE)
  say("stage C3 done")
} else {
  say("stage C3 skipped")
}

say("ALLDONE")
