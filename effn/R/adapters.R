# Input adapters. Four routes are supported and dispatched by class, with the
# Seurat route additionally dispatching on SeuratObject version.

#' Extract a counts matrix from a supported object
#'
#' @param object matrix, Matrix, data.frame, SingleCellExperiment or Seurat.
#' @param assay assay name (SingleCellExperiment / Seurat).
#' @param layer layer or slot name holding raw counts.
#' @return a genes x cells matrix-like object.
#' @keywords internal
#' @noRd
as_counts <- function(object, assay = NULL, layer = "counts") {
  if (is.data.frame(object)) object <- as.matrix(object)

  if (inherits(object, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Seurat object supplied but SeuratObject is not installed", call. = FALSE)
    }
    v <- utils::packageVersion("SeuratObject")
    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
    if (v >= "5.0.0") {
      out <- SeuratObject::GetAssayData(object, assay = assay, layer = layer)
    } else {
      out <- SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    }
    return(out)
  }

  if (inherits(object, "SingleCellExperiment")) {
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment supplied but SummarizedExperiment is not installed",
           call. = FALSE)
    }
    if (is.null(assay)) assay <- if ("counts" %in% SummarizedExperiment::assayNames(object)) "counts" else 1L
    return(SummarizedExperiment::assay(object, assay))
  }

  if (is.matrix(object) || inherits(object, "Matrix")) return(object)

  stop("unsupported object class: ", paste(class(object), collapse = "/"),
       "; supply a matrix, Matrix, data.frame, SingleCellExperiment or Seurat object",
       call. = FALSE)
}

#' Resolve a column name to the vector it names
#'
#' `group` is documented as a vector, but the first thing a user does with a
#' SingleCellExperiment or a Seurat object is pass the name of the column that
#' holds it. Accepting both costs one branch and removes a class of error
#' reports: `effn(sce, "group")` and `effn(sce, sce$group)` then agree.
#'
#' @param object the object passed to [effn()].
#' @param x a vector, or the name of one column of the sample metadata.
#' @return the vector.
#' @keywords internal
#' @noRd
resolve_column <- function(object, x) {
  if (!is.character(x) || length(x) != 1L) return(x)
  meta <- NULL
  if (inherits(object, "SingleCellExperiment") &&
      requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    meta <- as.data.frame(SummarizedExperiment::colData(object))
  } else if (inherits(object, "Seurat")) {
    meta <- object[[]]
  }
  if (!is.null(meta) && x %in% colnames(meta)) return(as.character(meta[[x]]))
  x
}

#' Guess donor labels from an object
#'
#' Four-tier route documented in DEFINITIONS.md section 4. Tier 1 reads donor
#' labels from sample metadata; tier 3 falls back to a sample-name heuristic and
#' flags the result; tier 4 stops.
#'
#' @param object the object passed to [effn()].
#' @param donor optional user-supplied donor vector.
#' @param sample optional sample-name vector used by the heuristic tier.
#' @return list with `donor`, `source`, `heuristic`.
#' @keywords internal
#' @noRd
extract_donor <- function(object, donor = NULL, sample = NULL) {
  if (!is.null(donor)) {
    return(list(donor = as.character(donor), source = "user", heuristic = FALSE))
  }

  meta <- NULL
  if (inherits(object, "SingleCellExperiment")) {
    meta <- as.data.frame(SummarizedExperiment::colData(object))
  } else if (inherits(object, "Seurat")) {
    meta <- object[[]]
  }

  if (!is.null(meta)) {
    ## "donor", "donor_id", "DonorID", "patient", "sample_patient" style names.
    ## The identifier suffix is optional and the whole token may be suffixed;
    ## matching only `<name>id` was a bug: the most common column of all, plain
    ## `donor`, fell through to the heuristic tier and was refused.
    cand <- grep("^(patient|donor|subject|individual|ind)([._]?id)?$",
                 colnames(meta), ignore.case = TRUE, value = TRUE)
    if (length(cand) > 0L) {
      return(list(donor = as.character(meta[[cand[1]]]),
                  source = paste0("metadata:", cand[1]), heuristic = FALSE))
    }
  }

  if (!is.null(sample)) {
    guess <- sub("^(.*?)[._-]?[0-9]+$", "\\1", as.character(sample))
    if (length(unique(guess)) == length(unique(sample)) && length(unique(guess)) > 1L) {
      return(list(donor = guess, source = "sample_name_heuristic", heuristic = TRUE))
    }
  }

  stop("donor labels unavailable: no donor field found in sample metadata and no ",
       "sample names supplied. Pass donor = <vector> explicitly; refusing to treat ",
       "the whole matrix as a single replicate.", call. = FALSE)
}
