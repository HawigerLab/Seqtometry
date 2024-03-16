#' Seqtometry scoring
#' @description
#' Scores signatures for the given matrix (which has already undergone QC, normalization, and imputation)
#' @param mat **matrix-like** Gene expression data (genes x cells)
#' @param signatures **named list of character** Signature genes (with same nomenclature system used in `mat`)
#' @param minmax **logical(1)** Whether to perform minmax transform on scoring results (default: TRUE)
#' @returns data.table (cells x signatures): single cell scores for each signature, where
#'   cell barcodes are stored in the "id" column
#' @importFrom data.table `:=`
#' @export
score <- function(mat, signatures, minmax = TRUE) {
  # Check parameter types
  checkmate::assert_multi_class(mat, c("matrix", "Matrix", "DelayedMatrix"))
  checkmate::assert_character(rownames(mat))
  checkmate::assert_list(signatures, type = "character")
  checkmate::assert_logical(minmax, len = 1L)

  # Means and standard deviations of all genes across all cells
  # for Z score transform performed inside Rcpp helper function
  mus <- MatrixGenerics::rowMeans2(mat, useNames = FALSE)
  sds <- MatrixGenerics::rowSds(mat, center = mus, useNames = FALSE)

  # Ignore any unexpressed genes
  nonzero <- which(sds != 0)
  if (length(nonzero) < nrow(mat)) {
    mat <- mat[nonzero, ]
    mus <- mus[nonzero]
    sds <- sds[nonzero]
  }

  # 0-based (for C++) row indices of signature genes
  gix <- .gene_indices(mat, signatures)

  # Loop (with opt-in parallelism) over all cells yielding signature score matrix
  # Calls Rcpp subroutine for weighted KS-like procedure
  # Convert to data.table for plotting later (in ggplot2)
  ret <- seq_len(ncol(mat)) |>
    future.apply::future_lapply(\(j) wks(mat[, j], gix, mus, sds)) |>
    data.table::transpose() |>
    data.table::setDT()

  # Apply minmax transform to scores
  if (minmax) ret[, names(ret) := lapply(ret, .minmax_scale)]

  # Assign cell barcodes and signature names to scores
  data.table::setnames(ret, names(signatures))
  if (!is.null(colnames(mat))) ret[, "id" := colnames(mat)]

  ret
}

#' Finds row indices of signature genes
#' @description For converting from character to integer based indexing
#' @param mat **matrix-like** Gene expression data (genes x cells)
#' @param gss **named list of character** Signature genes (with same nomenclature system as `mat`)
#' @returns **integer** 0-based indices (for passing to Rcpp function) of signature genes
#' @importFrom data.table `%chin%`
.gene_indices <- function(mat, gss) {
  n <- rownames(mat)
  lapply(gss, \(g) which(n %chin% g) - 1L) # Convert 1-based to 0-based indexing
}

#' Minmax transform
#' @description Scales input vector to unit range
#' @param x **numeric** Values to be scaled
#' @returns **numeric** Minmax transformed values
.minmax_scale <- function(x) {
  r <- range(x)
  (x - r[1]) / (r[2] - r[1])
}
