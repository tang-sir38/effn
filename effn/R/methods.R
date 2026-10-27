#' @export
print.effn_report <- function(x, ...) {
  s <- x$summary
  cat("<effn_report>\n")
  cat(sprintf("  design      : %s, %s cells, %d donors (%s)\n", s$design,
              format(s$nominal_n, big.mark = ","), s$n_donor,
              paste(names(s$n_donor_by_group), s$n_donor_by_group, sep = "=", collapse = ", ")))
  cat(sprintf("  donor source: %s%s\n", s$donor_source,
              if (isTRUE(s$heuristic)) "  [heuristic]" else ""))
  cat(sprintf(paste0("  genes       : %d tested, %d estimable, %d tech-dominated, ",
                     "%d low-resolution, %d low-expression, %d all-zero\n"),
              s$n_genes_tested, s$n_ok, s$n_unestimable[["tech_dominated"]],
              s$n_unestimable[["low_resolution"]],
              s$n_unestimable[["low_expression"]], s$n_unestimable[["all_zero"]]))
  cat(sprintf("  n_eff       : median %.2f of %d donors (95%% CI %.2f-%.2f)\n",
              s$median_n_eff, s$n_donor, s$median_n_eff_ci[1], s$median_n_eff_ci[2]))
  cat(sprintf("  MDE         : median %.2f log2 = %.2f-fold (95%% CI %.2f-%.2f log2)\n",
              s$median_mde_log2, s$median_mde_fold,
              s$median_mde_log2_ci[1], s$median_mde_log2_ci[2]))
  cat(sprintf("  loss        : cells->donors %.1fx, donors->n_eff %.2fx (depth share %.1f%%)\n",
              s$loss_donor, s$loss_depth, 100 * s$loss_depth_share))
  if (!is.null(s$n_required_2fold)) {
    two <- if (is.na(s$n_required_2fold)) "unreachable" else
      sprintf("%.0f (IQR %.0f-%.0f)", s$n_required_2fold,
              s$n_required_2fold_iqr[1], s$n_required_2fold_iqr[2])
    one5 <- if (is.na(s$n_required_1.5fold)) "unreachable" else
      sprintf("%.0f", s$n_required_1.5fold)
    cat(sprintf("  donors needed: %s for 2-fold, %s for 1.5-fold; %d present\n",
                two, one5, min(s$n_donor_by_group)))
  }
  cat(sprintf("  prescription: %s\n", s$prescription))
  invisible(x)
}

#' @export
summary.effn_report <- function(object, ...) {
  print(object)
  invisible(object$summary)
}

#' @export
as.data.frame.effn_report <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$per_gene
}
