#' Finds row indices of signature genes in given gene expression matrix
#' @param mat matrix-like: gene expression data (genes x cells)
#' @param gss named list of character: signature genes (with same nomenclature system as mat)
#' @returns 0-based indices (for passing to Rcpp function) of signature genes
#' @importFrom data.table `%chin%`
.gene_indices <- function(mat, gss) {
  n <- rownames(mat)
  lapply(gss, \(g) which(n %chin% g) - 1L) # Convert 1-based (from R) to 0-based indexing (for C++)
}

#' Scale vector to unit range = (0, 1)
#' @param x numeric: values to be scaled
#' @returns numeric: scaled values
.minmax_scale <- function(x) {
  r <- range(x)
  (x - r[1]) / (r[2] - r[1])
}

#' Seqtometry scoring
#'
#' Scores signatures for the given matrix (which has already undergone QC, normalization, and imputation)
#' @param mat matrix-like: gene expression data (genes x cells)
#' @param signatures named list of character: signature genes (with same nomenclature system used in mat)
#' @param minmax logical(1): whether to perform minmax transform on scoring results (default: TRUE)
#' @param fut_opt named list: additional options to customize the behavior of `%dofuture%` (default: NULL)
#' @returns matrix (cells x signatures): single cell scores for each signature
#' @importFrom MatrixGenerics rowMeans2 rowSds
#' @importFrom iterators icount
#' @importFrom doFuture `%dofuture%`
#' @importFrom foreach foreach
#' @importFrom data.table frankv
#' @details
#' Ensure that future globals are given enough memory.
#' For example, to grant a max of 16 GiB: `options(future.globals.maxSize = 16 * 1024^3)`
#' @export
score <- function(mat, signatures, minmax = TRUE, fut_opt = NULL) {
  gix <- .gene_indices(mat, signatures)      # 0-based (for C++) row indices of signature genes
  mus <- rowMeans2(mat)                      # Means of genes across all cells (for Z score transform)
  sds <- rowSds(mat, center = mus)           # Standard deviations of genes across all cells (for Z score transform)
  ctr <- ceiling(nrow(mat) / 2)              # Value for centering ranks
  aux <- c(ctr, ctr + (nrow(mat) %% 2 == 0)) # Centering value and start index for running sum

  # Loop (with opt-in parallelism) over all cells yielding signature score matrix
  ret <- foreach(i = icount(ncol(mat)), .combine = rbind, .options.future = fut_opt) %dofuture% {
    # Centered rank vector after Z score transform
    v <- frankv((mat[, i] - mus) / sds) - ctr

    # Loop over all signatures and perform weighted KS-like procedure (calls Rcpp function)
    sapply(gix, \(g) wks(v, g, aux))
  }

  # Apply minmax transform across all cells for each signature's scores
  if (minmax) apply(ret, 2, .minmax_scale) else ret
}
