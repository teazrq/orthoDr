#' @title Partial Sliced Averaged Variance Estimation
#' @name pSAVE
#' @export
#'
#' @description
#' Fits the partial SAVE model for dose-response settings. This model relies
#' on strong assumptions; its solution is primarily used as an initial value
#' in \code{\link{orthoDr_pdose}} optimization.
#'
#' @param x A finite numeric matrix (n x p) of features (continuous only),
#'   with at least two rows, one column, and full centred column rank.
#' @param a A finite numeric vector of observed dose levels, one per row of x.
#' @param r A finite numeric vector of rewards (outcomes), one per row of x.
#' @param ndr Integer between 1 and p; number of structural dimensions to estimate.
#'   Default: \code{2}.
#' @param nslices0 Positive integer; requested number of outcome slices for SAVE.
#'   Default: \code{2}.
#'
#' @return A matrix of unit-length directions in the original predictor
#'   coordinates, ordered by decreasing common-coordinate eigenvalue. Each
#'   column's largest absolute entry is positive (first entry in a tie).
#'   The columns need not be Euclidean-orthogonal.
#'
#' @details
#' Each central dose threshold contributes a grouped SAVE moment only when
#' every nonempty dose group has at least p + 1 observations and full centred
#' predictor rank. Repeated thresholds retain their repeated weights, and a
#' threshold with one observed dose group is allowed. Group and slice
#' covariances use their population divisors. Each group's symmetric-whitened
#' SAVE moment is transported to raw coefficient coordinates before pooling.
#' Its eigenproblem uses the group-size-weighted average within-group
#' covariance as a common reference, and the selected directions are mapped
#' back to the original predictor coordinates.
#'
#' Outcome slices of size one or two contribute only a scalar identity in
#' the common coordinates. They supply no directional information and are
#' excluded from the transported local moments without renormalizing the
#' remaining slice weights. This preserves their nondirectional fallback,
#' rather than transporting a local identity through unequal group covariances.
#' The construction corrects aggregation in incompatible local coordinates;
#' it does not reproduce historical pSAVE output. Affine predictor changes
#' transform the raw selected subspace accordingly when its eigenvalue
#' boundary is separated and numerical group admissibility is unchanged.
#' Statistical subspace recovery still requires appropriate SAVE assumptions.
#'
#' For p > 1, initialization fails if no usable dose threshold remains or if
#' every retained group has either one outcome slice or only slices smaller
#' than three. These cases supply no directional slice information. When used
#' to initialize \code{\link{orthoDr_pdose}}, supply \code{B.initial} to bypass
#' pSAVE if automatic initialization is unavailable. For a varying single
#' predictor and ndr = 1, the entire available feature subspace is returned
#' exactly as \code{matrix(1, 1, 1)} after input validation.
#'
#' @references
#' Feng, Z., Wen, X.M., Yu, Z. and Zhu, L. (2013). On partial sufficient
#' dimension reduction with applications to partially linear multi-index
#' models. \emph{Journal of the American Statistical Association},
#' 108(501), 237--246. \doi{10.1080/01621459.2012.746065}.
#'
#' @seealso \code{\link{orthoDr_pdose}}
#'
#' @examples
#' set.seed(1)
#' N <- 200; P <- 4
#' X <- matrix(rnorm(N * P), N, P)
#' dose <- runif(N, 0, 2)
#' reward <- X[, 1] + rnorm(N)
#' pSAVE(X, dose, reward, ndr = 1)

pSAVE <- function(x, a, r, ndr = 2, nslices0 = 2) {
  if (!is.matrix(x)) stop("X must be a matrix")
  if (!is.numeric(x) || is.complex(x)) stop("x must be numerical and real")
  if (nrow(x) != length(r) || nrow(x) != length(a)) stop("Number of observations do not match")
  if (nrow(x) < 2L || ncol(x) < 1L || any(!is.finite(x)))
    stop("x must have at least two rows and one column, with finite numeric values")
  for (name in c("a", "r")) {
    value <- get(name)
    if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
        any(!is.finite(value)))
      stop(name, " must be a finite numeric vector")
  }
  finite_scalar <- function(value) {
    is.numeric(value) && !is.complex(value) && length(value) == 1L &&
      is.null(dim(value)) && is.finite(value)
  }
  n <- nrow(x)
  p <- ncol(x)
  if (!finite_scalar(ndr) || ndr < 1 || ndr != floor(ndr) || ndr > p)
    stop("ndr must be an integer between 1 and ncol(x)")
  if (!finite_scalar(nslices0) || nslices0 < 1 || nslices0 != floor(nslices0))
    stop("nslices0 must be a positive integer")
  centred_rank <- function(z) qr(sweep(z, 2L, colMeans(z), "-"))$rank
  if (centred_rank(x) < p)
    stop("pSAVE requires full centred predictor rank; supply B.initial to orthoDr_pdose")
  if (p == 1L) return(matrix(1, 1, 1))

  # Preserve the central threshold multiset, including odd-size truncation and ties.
  Z <- sort(a)[seq.int(max(1L, floor(n / 2 - 50)), min(n, floor(n / 2 + 50)))]
  threshold_moments <- list()
  has_contrasts <- has_information <- FALSE
  for (cut in Z) {
    group <- as.numeric(a <= cut)
    indices <- lapply(unique(group), function(g) which(group == g))
    admissible <- vapply(indices, function(idx) {
      length(idx) >= p + 1L && centred_rank(x[idx, , drop = FALSE]) == p
    }, logical(1))
    if (!all(admissible)) next
    local <- .psave_grouped_moment(x, r, group, nslices0)
    threshold_moments[[length(threshold_moments) + 1L]] <- local
    has_contrasts <- has_contrasts || local$has.contrasts
    has_information <- has_information || local$has.information
  }
  if (!length(threshold_moments))
    stop("pSAVE has no admissible dose threshold; supply B.initial to orthoDr_pdose")
  if (!has_contrasts)
    stop("pSAVE has no outcome-slice contrasts; supply B.initial to orthoDr_pdose")
  if (!has_information)
    stop("pSAVE has no directional slice information: every retained group has one outcome slice or only slices smaller than three; supply B.initial to orthoDr_pdose")

  nt <- length(threshold_moments)
  reference <- raw_moment <- matrix(0, p, p)
  for (local in threshold_moments) {
    reference <- reference + local$covariance / nt
    raw_moment <- raw_moment + local$raw.moment / nt
  }
  roots <- .psave_covariance_roots(reference)
  common <- roots$sqrt %*% raw_moment %*% roots$sqrt
  # The omitted total fallback mass multiplies I and changes no directions.
  common <- (common + t(common)) / 2
  vectors <- eigen(common, symmetric = TRUE)$vectors[, seq_len(ndr), drop = FALSE]
  directions <- roots$inverse.sqrt %*% vectors
  directions <- sweep(directions, 2L, sqrt(colSums(directions^2)), "/")
  for (j in seq_len(ndr)) {
    pivot <- which.max(abs(directions[, j]))
    if (directions[pivot, j] < 0) directions[, j] <- -directions[, j]
  }
  directions
}

# Small explicit local moments avoid relying on the dr object's pooled layout.
.psave_covariance_roots <- function(covariance) {
  decomposition <- eigen((covariance + t(covariance)) / 2, symmetric = TRUE)
  if (any(!is.finite(decomposition$values)) || any(decomposition$values <= 0))
    stop("dose group covariance must be positive definite for whitening")
  vectors <- decomposition$vectors
  values <- sqrt(decomposition$values)
  list(sqrt = tcrossprod(sweep(vectors, 2L, values, "*"), vectors),
       inverse.sqrt = tcrossprod(sweep(vectors, 2L, values, "/"), vectors))
}

.psave_group_moment <- function(x, y, nslices) {
  n <- nrow(x)
  p <- ncol(x)
  centered <- sweep(x, 2L, colMeans(x), "-")
  covariance <- crossprod(centered) / n
  inverse <- .psave_covariance_roots(covariance)$inverse.sqrt
  z <- centered %*% inverse
  slices <- .psave_slices(y, nslices)
  moment <- matrix(0, p, p)
  fallback <- 0
  for (j in seq_len(slices$nslices)) {
    size <- slices$slice.sizes[j]
    if (size < 3L) {
      fallback <- fallback + size / n
      next
    }
    zj <- z[slices$slice.indicator == j, , drop = FALSE]
    zj <- sweep(zj, 2L, colMeans(zj), "-")
    difference <- diag(p) - crossprod(zj) / size
    moment <- moment + size / n * (difference %*% difference)
  }
  list(covariance = covariance, inverse.sqrt = inverse, moment = moment,
       fallback.mass = fallback, slices = slices)
}

.psave_grouped_moment <- function(x, y, group, nslices) {
  n <- nrow(x)
  p <- ncol(x)
  covariance <- raw_moment <- matrix(0, p, p)
  fallback <- 0
  has_contrasts <- has_information <- FALSE
  for (g in unique(group)) {
    selected <- group == g
    local <- .psave_group_moment(x[selected, , drop = FALSE], y[selected], nslices)
    weight <- sum(selected) / n
    covariance <- covariance + weight * local$covariance
    raw_moment <- raw_moment + weight *
      local$inverse.sqrt %*% local$moment %*% local$inverse.sqrt
    fallback <- fallback + weight * local$fallback.mass
    has_contrasts <- has_contrasts || local$slices$nslices > 1L
    has_information <- has_information ||
      (local$slices$nslices > 1L && any(local$slices$slice.sizes >= 3L))
  }
  list(covariance = covariance, raw.moment = raw_moment,
       fallback.mass = fallback, has.contrasts = has_contrasts,
       has.information = has_information)
}

# Repair only dr's confirmed one-breakpoint branch; preserve every other result.
.psave_slices <- function(y, h) {
  if (length(unique(y)) > h) {
    cty <- cumsum(table(y))
    first <- cty[which(cty >= floor(length(y) / h))[1L]]
    if (first >= length(y) - 2L)
      return(list(slice.indicator = rep(1, length(y)), nslices = 1L,
                  slice.sizes = length(y)))
  }
  dr::dr.slices(y, h)
}

