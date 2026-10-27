#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 13_quantile_gap.R
#
# How much the normal approximation understates the MDE at the donor counts
# single-cell studies actually use. The gap depends only on the design df, not
# on the data, so this is a table of quantiles rather than a re-analysis; it is
# kept as a file so the manuscript quotes a computed number.
#
# Writes 03_RESULTS/tables/t_normal_gap.tsv.
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "03_RESULTS/tables"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

p <- 0.80
z <- stats::qnorm(1 - 0.05 / 2) + stats::qnorm(p)

grid <- expand.grid(n_per_group = c(2L, 3L, 4L, 6L, 8L, 12L),
                    design = c("paired", "unpaired"),
                    stringsAsFactors = FALSE)
grid$df <- ifelse(grid$design == "paired",
                  grid$n_per_group - 1L,
                  grid$n_per_group + grid$n_per_group - 2L)
grid$t_sum <- stats::qt(1 - 0.05 / 2, grid$df) + stats::qt(p, grid$df)
grid$understatement_pct <- 100 * (1 - z / grid$t_sum)
grid <- grid[order(grid$design, grid$n_per_group), ]

write.table(grid, file.path(outdir, "t_normal_gap.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
print(grid[, c("design", "n_per_group", "df", "t_sum", "understatement_pct")],
      row.names = FALSE, digits = 4)
cat("\nwrote", file.path(outdir, "t_normal_gap.tsv"), "\n")
cat("QUANTILEGAPDONE\n")
