#' Performs MAGIC imputation
#' @description
#' Calculates a graph diffusion operator for the given input matrix and applies it to produce an imputed matrix.
#' @param gex **matrix or Matrix** Gene expression values (that has passed quality control).
#' @param transpose **logical(1)** Whether to transpose gex (make it cells x genes) prior to downstream operations.
#' @param do_norm **logical(1)** Whether to perform LogCP10K normalization on gex.
#' @param pca **matrix (cells x PCs) or NULL** Precomputed principal component matrix (or NULL to derive it from gex).
#' @param npc **integer(1)** Number of principal components (min = 1) to calculate.
#' @param scale **logical(1)** Whether to scale columns of input matrix (after centering) to unit variance prior to PCA.
#' @param knn **integer(1)** Number of nearest neighbors (min = 2) to consider during distance calculation.
#' @param ka **integer(1)** Number of nearest neighbors (min = 2, max <= knn) to use for the adaptive kernel.
#' @param dist_metric **character(1)** Type of metric to use for distance calculations during kNN search.
#' @param dft **NULL or integer(1)** Automatic (NULL) or user-defined (integer) diffusion time (min = 1, max = 16).
#' @param t_max **integer(1)** Maximum diffusion time to test when using automatic diffusion time (min = 1, max = 16).
#' @param tol **numeric(1)** Threshold for Procrustes disparity (min = 0, max = 1) between successive diffusion times.
#' @param exact_solver **logical(1)** Whether to perform imputation in gene space (TRUE) or PCA space (FALSE).
#' @param conserve_memory **logical(1)** Whether to avoid allocating a large dense matrix when `exact_solver = FALSE`.
#' @param env_ret **logical(1)** Return all variables in the environment (TRUE) or just the imputed matrix (FALSE).
#' @param verbose **logical(1)** Whether to issue print statements at different major parts of the algorithm.
#' @returns **matrix-like or list**
#' If `env_ret = FALSE`, then just the imputed matrix.
#' Otherwise the function environment as a list containing all parameters (possibly modified) as well as
#' - imp **matrix or DelayedMatrix** Imputed matrix.
#' - aff **dgCMatrix** Markov affinity matrix (graph diffusion operator).
#' - pca **list** Possibly computed (if pca was NULL), yielding a four element list, where:
#'     - `x` **matrix (cells x PCs)** The principal components matrix (scaled left singular vectors).
#'     - `v` **matrix (genes x PCs)** The rotation matrix (right singular vectors).
#'     - `center` **integer (cells)** The centering vector.
#'     - `scale` **integer (cells) or NULL** The scaling vector (or NULL if no scaling was applied).
#' @importFrom zeallot `%<-%`
#' @export
impute <- function(gex,
    transpose = TRUE, do_norm = FALSE,
    pca = NULL, npc = 100L, scale = TRUE,
    knn = 16L, ka = 6L, dist_metric = "euclidean",
    dft = NULL, t_max = 16L, tol = 0.001,
    exact_solver = TRUE, conserve_memory = FALSE,
    env_ret = FALSE, verbose = TRUE) {
  .check_params(environment())

  if (transpose) {
    if (verbose) cat("Transposing input matrix\n")
    gex <- Matrix::t(gex)
  }
  if (do_norm) {
    if (verbose) cat("Normalizing input matrix\n")
    gex <- .normalize(gex)
  }

  pca <- if (is.null(pca)) {
    if (verbose) cat("Performing PCA\n")
    .calc_pca(gex, npc, scale)
  } else {
    if (verbose) cat("Skipping PCA (custom PC matrix provided)\n")
    list(x = pca, v = NULL, center = NULL, scale = NULL)
  }

  if (verbose) cat("Calculating diffusion operator\n")
  aff <- .calc_diff_op(pca$x, knn, ka, dist_metric)

  if (verbose) cat("Applying diffusion operator")
  c(imp, dft) %<-% .apply_diff_op(gex, pca$x, aff, dft, t_max, tol, exact_solver)
  imp <- `if`(exact_solver,
    Matrix::t(imp) |> as.matrix(),
    .invert_pca(imp, pca$v, pca$center, pca$scale, conserve_memory))
  dimnames(imp) <- dimnames(gex) |> rev()

  if (env_ret) as.list(environment()) else imp
}

#' Parameter validation
#' @description Checks that all the parameters used in MAGIC imputation are permissible.
#' @param args The arguments to the magic_impute function
#' @returns NULL (but stops execution for invalid parameters)
.check_params <- function(args) {
  checkmate::assert_multi_class(args$gex, c("Matrix", "matrix"))
  checkmate::assert_logical(args$transpose, len = 1L)
  checkmate::assert_logical(args$do_norm, len = 1L)
  checkmate::assert_logical(args$env_ret, len = 1L)
  checkmate::assert_logical(args$verbose, len = 1L)
  checkmate::assert_logical(args$exact_solver, len = 1L)
  checkmate::assert_logical(args$conserve_memory, len = 1L)
  if (is.null(args$pca)) {
    checkmate::assert_int(args$npc, lower = 1L, upper = nrow(args$gex))
    checkmate::assert_logical(args$scale, len = 1L)
  } else {
    checkmate::assert_true(args$exact_solver)
    checkmate::assert_matrix(args$pca)
    checkmate::assert_true(nrow(args$pca) == `if`(args$transpose, ncol, nrow)(args$gex))
  }
  checkmate::assert_int(args$knn, lower = 2L)
  checkmate::assert_int(args$ka, lower = 2L, upper = args$knn)
  checkmate::assert_true(args$dist_metric %in% c("euclidean", "l2", "cosine", "ip"))
  if (is.null(args$dft)) {
    checkmate::assert_int(args$t_max, lower = 1L, upper = 16L)
    checkmate::assert_number(args$tol, lower = 0, upper = 1)
  } else {
    checkmate::assert_int(args$dft, lower = 1L, upper = 16L)
  }
  NULL
}

#' LogCP10K transform
#' @description Simple normalization method for scRNA-seq data
#' @param gex Gene expression matrix (cells x genes)
#' @returns **matrix or dgCMatrix** Transformed (normalized) matrix
.normalize <- function(gex) {
  scale_factors <- 1e4 / Matrix::rowSums(gex)
  log1p(gex * scale_factors)
}

#' PCA wrapper
#' @description Calculate leading principal components (via truncated singular value decomposition).
#' @param gex **matrix or Matrix** Gene expression matrix
#' @param npc **numeric(1)** Number of leading principal components to compute
#' @param scale **logical(1)** Whether to scale genes to unit variance
#' @returns **list** PC loading/rotation matrices as well as centering/scaling vectors.
.calc_pca <- function(gex, npc, scale) {
  ctr <- MatrixGenerics::colMeans2(gex, useNames = FALSE)
  sdv <- if (scale) MatrixGenerics::colSds(gex, center = ctr, useNames = FALSE)
  svd <- RSpectra::svds(gex, npc, opts = list(center = ctr, scale = sdv))
  sweep(svd$u, 2, svd$d * sqrt(nrow(gex) - 1), "*") |>
    list(x = _, v = svd$v, center = ctr, scale = sdv)
}

#' Compute diffusion operator
#' @description Calculate graph diffusion operator (Markov affinity matrix).
#' @param pcs **matrix** Principal components matrix (used for kNN search)
#' @param knn **integer(1)** Number of nearest neighbors to search for
#' @param ka **integer(1)** Number of nearest neighbors to use for adaptive kernel
#' @param dist_metric **character(1)** Type of metric to use for distance calculations during kNN search.
#' @returns **dgCMatrix** Markov affinity matrix.
.calc_diff_op <- function(pcs, knn, ka, dist_metric) {
  nbr <- RcppHNSW::hnsw_knn(pcs, knn, distance = dist_metric)
  # Distances to affinities using Gaussian kernel with adaptive width
  aff <- Matrix::sparseMatrix(
    i = nrow(pcs) |> seq_len() |> rep(knn), j = as.vector(nbr$idx),
    x = exp(-(as.vector(nbr$dist) / nbr$dist[, ka])^2),
    dims = c(nrow(pcs), nrow(pcs)), dimnames = list(rownames(pcs), rownames(pcs)))
  aff <- aff + Matrix::t(aff) # Symmetrization
  aff / Matrix::rowSums(aff)  # Markov normalization
}

#' Perform data diffusion
#' @description Apply diffusion operator to data in order to perform imputation.
#' @param gex **matrix or Matrix** Gene expression matrix
#' @param pcs **matrix** Principal components matrix
#' @param aff **dgCMatrix** Markov affinity matrix
#' @param dft **NULL or integer(1)** Diffusion time
#' @param t_max **integer(1)** Maximum diffusion time
#' @param tol **numeric(1)** Tolerance for Procrustes disparity
#' @param exact_solver **logical(1)** Perform imputation in gene space
#' @returns **list(matrix, integer(1))** Imputed matrix and diffusion time used.
.apply_diff_op <- function(gex, pcs, aff, dft, t_max, tol, exact_solver) {
  # (aff %^% pow) %*% x; use repeated sparse matrix multiplication instead of exponentiation
  mtx_mul <- \(x, pow) purrr::reduce(seq_len(pow), \(m, i) aff %*% m, .init = x)
  imp <- if (!is.null(dft)) {
    `if`(exact_solver, gex, pcs) |> mtx_mul(dft)
  } else { # Determine optimal diffusion time using the PC matrix (for performance reasons)
    m <- pcs
    stop <- FALSE
    dft <- 0L
    while (!stop && dft <= t_max) {
      m0 <- m
      m <- aff %*% m0
      stop <- .procrustes(m0, m) < tol
      dft <- dft + 1L
    }
    `if`(exact_solver, mtx_mul(gex, dft), m)
  }
  list(imp, dft)
}

#' Procrustes disparity
#' @description Calculates symmetric Procrustes distance (adapted from MATLAB procrustes).
#' @param x **matrix**
#' @param y **matrix**
#' @returns **numeric(1)** Procrustes disparity between input matrices.
.procrustes <- function(x, y) {
  pre <- \(m) {
    m <- scale(m, scale = FALSE)
    m / norm(m, "F")
  }
  tr <- crossprod(pre(x), pre(y)) |> svd() |> _$d |> sum()
  1 - tr^2
}

#' PCA inversion
#' @description Reverses operations done for PCA: back-rotation, unscaling, and uncentering.
#' @param pcs **matrix** The principal components (scaled left singular vectors).
#' @param rot **matrix** The rotation matrix (right singular vectors).
#' @param ctr **integer** The centering vector.
#' @param sdv **integer or NULL** The scaling vector (or NULL if no scaling was applied).
#' @param low_mem **logical(1)** Whether to use delayed operations to reduce memory usage.
#' @returns **matrix or DelayedMatrix** `rot %*% t(pcs) * sdv + ctr`
.invert_pca <- function(pcs, rot, ctr, sdv, low_mem) {
  ret <- `if`(low_mem, BiocSingular::LowRankMatrix, tcrossprod)(rot, pcs)
  if (!is.null(sdv)) ret <- ret * sdv
  ret + ctr
}