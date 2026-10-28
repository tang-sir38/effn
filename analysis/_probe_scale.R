## Helper, not part of the numbered pipeline: run after 03_depth_and_export.R.
## Reads depth_core_genes.txt and gse96583_singlets.rds from $EFFN_OUT.
##
## Same depth curve on two scales.
## Var(observed donor-level value) = Var(biological) + Var(technical) is an exact
## additive decomposition only on a linear scale. On log2(CPM + 1) the transform
## bends the distribution by an amount that depends on the count rate, so the
## residual sigma_bio inherits a depth trend. This script measures both.
suppressPackageStartupMessages({ library(Matrix); library(effn) })

OUT <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
ob <- readRDS(file.path(OUT, "gse96583_singlets.rds"))
M <- as(methods::as(ob$counts, "CsparseMatrix"), "dgCMatrix")
donor <- ob$donor; grp <- ob$group
core <- readLines(file.path(OUT, "depth_core_genes.txt"))
M <- M[core, , drop = FALSE]
cat("core genes:", nrow(M), "\n")

pb_cpm <- function(counts, key) {
  ug <- unique(key)
  out <- matrix(0, length(ug), nrow(counts), dimnames = list(ug, rownames(counts)))
  for (d in ug) out[d, ] <- as.numeric(Matrix::rowSums(counts[, which(key == d), drop = FALSE]))
  sweep(out, 1, rowSums(out), "/") * 1e6
}

tech_boot <- function(counts, key, B = 40L, log = TRUE, seed = 20261021L) {
  set.seed(seed)
  ug <- unique(key); by <- split(seq_len(ncol(counts)), key)
  tm <- matrix(NA_real_, length(ug), nrow(counts), dimnames = list(ug, rownames(counts)))
  for (d in ug) {
    ci <- by[[d]]; Nd <- length(ci); if (Nd < 4L) next
    w <- rmultinom(B, Nd, rep(1 / Nd, Nd))
    W <- sparseMatrix(i = rep(seq_len(Nd), times = B), j = rep(seq_len(B), each = Nd),
                      x = as.numeric(w), dims = c(Nd, B))
    s <- as.matrix(counts[, ci, drop = FALSE] %*% W)
    tot <- colSums(s)
    y <- if (log) log2(sweep(s, 2, tot, "/") * 1e6 + 1) else sweep(s, 2, tot, "/") * 1e6
    mu <- rowMeans(y); tm[d, ] <- rowSums((y - mu)^2) / (B - 1L)
  }
  tm
}

rows <- list()
for (cc in c(10, 25, 50, 100, 250, 500, 1000)) {
  set.seed(20261021L)
  keep <- unlist(lapply(split(seq_along(donor), donor), function(ii)
    sample(ii, min(cc, length(ii)))), use.names = FALSE)
  Ms <- M[, keep, drop = FALSE]; dd <- donor[keep]; gg <- grp[keep]
  key <- paste(dd, gg, sep = "\r")
  for (sc in c("log", "linear")) {
    logsc <- identical(sc, "log")
    y <- pb_cpm(Ms, key)
    if (logsc) y <- log2(y + 1)
    kd <- sub("\r.*$", "", rownames(y)); kg <- sub("^.*\r", "", rownames(y))
    ord <- order(kd, kg); yw <- y[ord, , drop = FALSE]
    dif <- yw[kg[ord] == "stim", , drop = FALSE] - yw[kg[ord] == "ctrl", , drop = FALSE]
    s2obs <- apply(dif, 2, var)
    tm <- tech_boot(Ms, key, B = 40L, log = logsc)
    t1 <- tm[kg == "ctrl", , drop = FALSE]; t2 <- tm[kg == "stim", , drop = FALSE]
    rownames(t1) <- kd[kg == "ctrl"]; rownames(t2) <- kd[kg == "stim"]
    s2tech <- colMeans(t1[rownames(t2), , drop = FALSE] + t2, na.rm = TRUE)
    s2bio <- s2obs - s2tech
    rho <- ifelse(s2bio > 0, s2tech / s2bio, NA_real_)
    ok <- is.finite(rho) & rho >= 0
    rows[[length(rows) + 1L]] <- data.frame(
      cells = cc, scale = sc, n_ok = sum(ok),
      med_sigma_tech = median(sqrt(s2tech[ok])), med_sigma_bio = median(sqrt(s2bio[ok])),
      med_rho = median(rho[ok]), med_share = median(rho[ok] / (1 + rho[ok])),
      med_n_eff = median(8 / (1 + rho[ok]))
    )
    cat(sprintf("cells %4d %-6s ok %3d sigma_tech %.3f sigma_bio %.3f rho %.3f share %.3f n_eff %.2f\n",
                cc, sc, sum(ok), median(sqrt(s2tech[ok])), median(sqrt(s2bio[ok])),
                median(rho[ok]), median(rho[ok] / (1 + rho[ok])), median(8 / (1 + rho[ok]))))
  }
}
out <- do.call(rbind, rows)
utils::write.table(out, file.path(OUT, "depth_curve_two_scales.tsv"), sep = "\t",
                   quote = FALSE, row.names = FALSE)
cat("SCALEDONE\n")
