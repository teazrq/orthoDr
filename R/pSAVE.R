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
#' @return A matrix whose columns are the estimated basis vectors of the
#'   central subspace, ordered by eigenvalues.
#'
#' @details
#' Each central dose threshold contributes a grouped SAVE moment only when
#' every nonempty dose group has at least p + 1 observations and full centred
#' predictor rank. Repeated thresholds retain their repeated weights, and a
#' threshold with one observed dose group is allowed. Outcome slices of size
#' one or two retain the SAVE dependency's identity contribution.
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
#' Feng, Z., Wen, M.X., Yu, Z. and Zhu, L. (2013). On partial sufficient
#' dimension reduction with applications to partially linear multi-index
#' models. \emph{Journal of the American Statistical Association},
#' 108(501), 237--256. \doi{10.1080/01621459.2013.849167}.
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
  groups <- list()
  has_contrasts <- has_information <- FALSE
  for (cut in Z) {
    group <- as.numeric(a <= cut)
    indices <- lapply(unique(group), function(g) which(group == g))
    admissible <- vapply(indices, function(idx) {
      length(idx) >= p + 1L && centred_rank(x[idx, , drop = FALSE]) == p
    }, logical(1))
    if (!all(admissible)) next
    groups[[length(groups) + 1L]] <- group
    slices <- lapply(indices, function(idx) .psave_slices(r[idx], nslices0))
    has_contrasts <- has_contrasts || any(vapply(slices, function(s) {
      s$nslices > 1L
    }, logical(1)))
    has_information <- has_information || any(vapply(slices, function(s) {
      s$nslices > 1L && any(s$slice.sizes >= 3L)
    }, logical(1)))
  }
  if (!length(groups))
    stop("pSAVE has no admissible dose threshold; supply B.initial to orthoDr_pdose")
  if (!has_contrasts)
    stop("pSAVE has no outcome-slice contrasts; supply B.initial to orthoDr_pdose")
  if (!has_information)
    stop("pSAVE has no directional slice information: every retained group has one outcome slice or only slices smaller than three; supply B.initial to orthoDr_pdose")

  train <- list(x = x, a = a, r = r)
  M_total <- matrix(0, p, p)
  for (group in groups) {
    train$a <- group
    dimdr <- dr(formula = r ~ x, data = train, group = ~a,
                nslices = nslices0, chi2approx = "bx", numdir = p,
                method = "save", slice.function = .psave_slices)
    M_total <- M_total + dimdr$M
  }
  Beta <- eigen(M_total, symmetric = FALSE, only.values = FALSE, EISPACK = FALSE)
  Beta$vectors[, seq_len(ndr), drop = FALSE]
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

