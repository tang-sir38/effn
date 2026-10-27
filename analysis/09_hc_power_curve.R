#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# 09_hc_power_curve.R  (cloud host)
# Same simulated data summaries, larger designs: does Hierarchicell's hurdle
# power rise with donors the way the normal-theory calculation says it should?
# Produces out/hc_power_curve.tsv for the Supplementary comparison.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(hierarchicell)
})
OUT <- Sys.getenv("EFFN_OUT", file.path(Sys.getenv("EFFN_ROOT", "/root/effn"), "out"))
ds <- readRDS(file.path(OUT, "hc_data_summaries.rds"))
say <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")
  flush(stdout())
}

grid <- c(8L, 16L, 32L)
rows <- list()
for (n in grid) {
  msgs <- character(0)
  grab <- function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }
  t0 <- Sys.time()
  set.seed(20261021)
  invisible(withCallingHandlers(
    power_hierarchicell(ds, n_genes = 100, n_per_group = n,
                        cells_per_case = 300, cells_per_control = 300,
                        pval = 0.05, foldchange = 2),
    message = grab))
  el <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2)
  pick <- function(lab) {
    v <- grep(lab, msgs, value = TRUE)
    if (!length(v)) return(NA_real_)
    as.numeric(sub("^.*: ", "", v[length(v)]))
  }
  rows[[length(rows) + 1L]] <- data.frame(
    n_per_group = n, n_genes = 100L, cells_per_donor = 300L, foldchange = 2,
    hurdle_power = pick("Hurdle Power"),
    continuous_power = pick("Continuous Power"),
    discrete_power = pick("Discrete Power"),
    minutes = el)
  say("n=", n, " hurdle=", rows[[length(rows)]]$hurdle_power, " (", el, " min)")
}
res <- do.call(rbind, rows)
write.table(res, file.path(OUT, "hc_power_curve.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
print(res)
say("HCCURVEDONE")
