#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 01_build_gse96583.R  (run on the cloud host)
# Kang 2018 GSE96583 batch 2: build the count matrix, attach donor labels from
# the authors' per-cell demuxlet table, run effn(), and export everything the
# local pipeline needs.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(effn)
})

ROOT <- file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "GSE96583")
EX   <- file.path(ROOT, "extract")
OUT  <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ---- 1. genes ------------------------------------------------------------
genes <- utils::read.delim(file.path(ROOT, "GSE96583_batch2.genes.tsv.gz"),
                           header = FALSE, stringsAsFactors = FALSE)
cat("genes:", nrow(genes), "\n")

read_bc <- function(f) utils::read.delim(f, header = FALSE, stringsAsFactors = FALSE)[[1]]

## ---- 2. the two batch-2 libraries ---------------------------------------
ctrl_m <- Matrix::readMM(file.path(EX, "GSM2560248_2.1.mtx.gz"))
stim_m <- Matrix::readMM(file.path(EX, "GSM2560249_2.2.mtx.gz"))
ctrl_b <- read_bc(file.path(EX, "GSM2560248_barcodes.tsv.gz"))
stim_b <- read_bc(file.path(EX, "GSM2560249_barcodes.tsv.gz"))
cat("ctrl:", nrow(ctrl_m), "x", ncol(ctrl_m), " stim:", nrow(stim_m), "x", ncol(stim_m), "\n")
stopifnot(nrow(ctrl_m) == nrow(genes), nrow(stim_m) == nrow(genes),
          length(ctrl_b) == ncol(ctrl_m), length(stim_b) == ncol(stim_m))

## 313 barcodes occur in both libraries. The authors disambiguated them in the
## demuxlet table by appending "-1" to the stimulated library (suffix becomes
## "-11"). Verified: all 313 such rows carry stim in the meta table.
ov <- intersect(stim_b, ctrl_b)
cat("barcodes in both libraries:", length(ov), "\n")
stim_b2 <- ifelse(stim_b %in% ov, paste0(stim_b, "1"), stim_b)

M <- cbind(ctrl_m, stim_m)
bc <- c(ctrl_b, stim_b2)
rownames(M) <- genes[[1]]
colnames(M) <- bc
cat("combined:", nrow(M), "x", ncol(M), " duplicated:", anyDuplicated(bc), "\n")
stopifnot(anyDuplicated(bc) == 0)

## ---- 3. donor labels -----------------------------------------------------
## the header line has 7 fields while the data rows have 8 (the barcode is an
## unnamed first column), so read it manually instead of trusting read.delim
raw <- readLines(gzfile(file.path(ROOT, "b2_tsne.tsv.gz")))
hdr <- strsplit(raw[1], "\t")[[1]]
meta <- utils::read.delim(text = paste(raw[-1], collapse = "\n"), header = FALSE,
                          stringsAsFactors = FALSE)
names(meta) <- c("barcode", hdr)
cat("meta:", nrow(meta), "cols:", paste(names(meta), collapse = ","), "\n")
idx <- match(bc, meta$barcode)
cat("matched:", sum(!is.na(idx)), "of", length(bc), "\n")
stopifnot(!anyNA(idx))          # one-to-one: every cell has a row
donor <- meta$ind[idx]
stim  <- meta$stim[idx]
mult  <- meta$multiplets[idx]
ct    <- meta$cell[idx]

sing <- mult == "singlet"
cat("singlets:", sum(sing), "of", length(sing), "\n")
M <- M[, sing]; bc <- bc[sing]; donor <- donor[sing]
stim <- stim[sing]; ct <- ct[sing]

cat("final:", nrow(M), "genes x", ncol(M), "cells\n")
print(table(donor, stim))
stopifnot(all(table(donor, stim) > 0))

saveRDS(list(counts = M, donor = donor, group = stim, celltype = ct, barcode = bc),
        file.path(OUT, "gse96583_singlets.rds"), compress = "xz")
utils::write.table(data.frame(barcode = bc, donor = donor, group = stim, celltype = ct),
                   gzfile(file.path(OUT, "gse96583_cells.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 4. effn -------------------------------------------------------------
cat("running effn ...\n")
t0 <- Sys.time()
rep1 <- effn(M, group = stim, donor = donor, n_boot = 200L, n_split = 20L,
             seed = 20261021L)
cat("effn took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
print(rep1)

utils::write.table(rep1$per_gene, gzfile(file.path(OUT, "gse96583_per_gene.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)
s <- rep1$summary
writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, digits = 10, pretty = TRUE),
           file.path(OUT, "gse96583_summary.json"))

## ---- 5. determinism ------------------------------------------------------
rep2 <- effn(M, group = stim, donor = donor, n_boot = 200L, n_split = 20L,
             seed = 20261021L)
same <- identical(rep1$per_gene, rep2$per_gene) && identical(rep1$summary, rep2$summary)
cat("DETERMINISM identical:", same, "\n")
writeLines(paste0("determinism_identical=", same), file.path(OUT, "determinism.txt"))

## ---- 6. depth curve by empirical downsampling ---------------------------
set.seed(20261021)
min_cells <- min(table(donor))
ds <- sort(unique(pmin(c(10, 25, 50, 100, 250, 500, 1000, 2500), min_cells)))
rows <- list()
for (cc in ds) {
  keep <- unlist(lapply(split(seq_along(donor), donor), function(ii)
    sample(ii, min(cc, length(ii)))))
  Ms <- M[, keep]; dd <- donor[keep]; gg <- stim[keep]
  y <- effn:::pb_norm(Ms, paste(dd, gg, sep = "\r"))
  kd <- sub("\r.*$", "", rownames(y)); kg <- sub("^.*\r", "", rownames(y))
  ord <- order(kd, kg)
  yw <- y[ord, , drop = FALSE]
  dif <- yw[kg[ord] == "stim", , drop = FALSE] - yw[kg[ord] == "ctrl", , drop = FALSE]
  s2 <- apply(dif, 2, stats::var)
  tc <- effn:::sigma_tech_splithalf(Ms, paste(dd, gg, sep = "\r"), K = 10L, seed = 20261021L)
  tm <- tc$t_dg
  tk <- sub("\r.*$", "", rownames(tm)); tg <- sub("^.*\r", "", rownames(tm))
  t1 <- tm[tg == "ctrl", , drop = FALSE]; t2 <- tm[tg == "stim", , drop = FALSE]
  rownames(t1) <- tk[tg == "ctrl"]; rownames(t2) <- tk[tg == "stim"]
  t2 <- t2[rownames(t1), , drop = FALSE]
  ttot <- t1 + t2
  tech <- colMeans(ttot, na.rm = TRUE)
  bio <- s2 - tech
  rho <- ifelse(bio > 0, tech / bio, NA_real_)
  rows[[length(rows) + 1]] <- data.frame(
    cells_per_donor = cc,
    median_rho = median(rho, na.rm = TRUE),
    depth_share = median(rho / (1 + rho), na.rm = TRUE),
    frac_tech_dominated = mean(!is.finite(rho) | rho <= 0),
    n_cells_used = length(keep)
  )
  cat("depth", cc, "share", round(rows[[length(rows)]]$depth_share, 4), "\n")
}
dc <- do.call(rbind, rows)
utils::write.table(dc, file.path(OUT, "depth_curve.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
print(dc)

## ---- 7. subsample for local re-checks -----------------------------------
set.seed(20261022)
sub <- unlist(lapply(split(seq_along(donor), donor), function(ii)
  sample(ii, min(100, length(ii)))))
Ms <- M[, sub]
tr <- summary(Ms)
utils::write.table(data.frame(i = tr$i, j = tr$j, x = tr$x),
                   gzfile(file.path(OUT, "subsample_triplets.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
writeLines(rownames(Ms), gzfile(file.path(OUT, "subsample_genes.txt.gz")))
utils::write.table(data.frame(barcode = colnames(Ms), donor = donor[sub], group = stim[sub]),
                   gzfile(file.path(OUT, "subsample_cells.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 8. pseudobulk for figure checks ------------------------------------
pb <- effn:::pb_norm(M, paste(donor, stim, sep = "|"))
utils::write.table(as.data.frame(pb), gzfile(file.path(OUT, "pseudobulk_log2cpm.tsv.gz")),
                   sep = "\t", quote = FALSE, row.names = TRUE)

cat("ALLDONE\n")
