#' Personalized Dose Estimation via Dimension Reduction
#'
#' @description
#' Fits a personalized dose model using the direct learning or pseudo-direct
#' learning method of Zhou, Zhu & Zeng (2021). Direct learning targets a dose
#' rule, whereas pseudo-direct learning first targets the conditional mean
#' reward surface. The direction matrix \code{B} therefore has a different
#' target under the two methods.
#'
#' @name orthoDr_pdose
#' @export
#'
#' @param x A numeric matrix of features (continuous only).
#' @param a A numeric vector of observed dose levels.
#' @param r A numeric vector of observed rewards (outcomes). Larger values
#'   are preferred.
#' @param ndr Number of directions to estimate.
#' @param B.initial Initial \code{B} matrix. If \code{NULL} (default),
#'   the partial SAVE model (\code{\link{pSAVE}}) is used. If specified,
#'   must have \code{nrow(B.initial) == ncol(x)} and
#'   \code{ncol(B.initial) == ndr}. Will be orthogonalized via Gram-Schmidt.
#' @param bw Kernel bandwidth, assuming each variable has unit variance.
#'   If \code{NULL} (default), Silverman's rule is used.
#' @param lambda Penalty level for kernel ridge regression. If a vector of
#'   values is provided, GCV is used to select the best value.
#'   Default: \code{0.1}.
#' @param K Number of candidate dose grid points used during fitting. The grid
#'   spans the observed dose range. Must be a finite scalar at least one;
#'   fractional values are rounded down once. Default: \code{sqrt(length(r))}.
#'   A single grid point uses the minimum observed dose. Predictions are
#'   continuous and are not constrained to the grid or observed dose range.
#' @param method The method: \code{"direct"} or \code{"pseudo_direct"}.
#'   Default: \code{"direct"}.
#' @param keep.data Logical; if \code{TRUE}, the original data is stored in
#'   the returned object for prediction. Default: \code{FALSE}.
#' @param control A list of optimizer tuning parameters (see
#'   \code{\link{ortho_optim}} for details): \code{rho}, \code{eta},
#'   \code{gamma}, \code{tau}, \code{epsilon}, \code{btol}, \code{ftol},
#'   \code{gtol}.
#' @param maxitr Maximum number of iterations. Default: \code{500}.
#' @param verbose Integer; if > 0, prints iteration progress. Default: 0.
#' @param ncore Integer; number of cores for parallel computation via OpenMP.
#'   0 = auto-detect. Default: \code{0}.
#'
#' @return An object of class \code{c("orthoDr", "fit", "pdose")}, which is a
#'   list containing:
#'   \describe{
#'     \item{\code{B}}{The estimated direction matrix (columns are orthonormal).}
#'     \item{\code{fn}}{The final objective function value.}
#'     \item{\code{itr}}{Number of iterations performed.}
#'     \item{\code{converge}}{Logical; \code{TRUE} when a convergence criterion
#'       was reached. For the direct method, the gradient criterion is
#'       conditional on the current policy coefficients.}
#'     \item{\code{method}}{The method used (\code{"direct"} or
#'       \code{"pseudo_direct"}).}
#'     \item{\code{keep.data}}{Whether original data was retained.}
#'   }
#'
#' @details
#' With \code{method = "direct"}, \code{B} is chosen to improve estimated mean
#' reward under a dose rule based on \eqn{B^\top x}. With
#' \code{method = "pseudo_direct"}, \code{B} is chosen by fitting the
#' conditional mean reward as a function of dose and \eqn{B^\top x}, then a
#' dose rule is constructed. The methods need not estimate the same subspace.
#'
#' The dose grid is used to select candidate doses during fitting. Both
#' methods fit a kernel ridge regression rule to the selected doses.
#' Prediction evaluates this continuous rule, so predicted doses need not lie
#' on the grid or within the observed dose range.
#'
#' Dose comparisons use the standard deviation of the observed training doses
#' for both observed and candidate doses. Returned doses remain in their
#' original units. The direct method updates policy coefficients during line
#' search but computes its gradient with those coefficients held fixed; its
#' stopping criterion does not certify stationarity of a profiled objective.
#'
#' @references
#' Zhou, W., Zhu, R. and Zeng, D. (2021). A parsimonious personalized dose-finding model via
#' dimension reduction. \emph{Biometrika}, 108(3), 643--659.
#' \doi{10.1093/biomet/asaa087}.
#'
#' @examples
#' # generate personalized dose data
#' exampleset <- function(size, ncov) {
#'   X <- matrix(runif(size * ncov, -1, 1), ncol = ncov)
#'   A <- runif(size, 0, 2)
#'   Edr <- as.matrix(c(0.5, -0.5))
#'   D_opt <- X %*% Edr + 1
#'   mu <- 2 + 0.5 * (X %*% Edr) - 7 * abs(D_opt - A)
#'   R <- rnorm(length(mu), mu, 1)
#'   R <- R - min(R)
#'   list(X = X, A = A, R = R, D_opt = D_opt, mu = mu)
#' }
#'
#' set.seed(123)
#' n <- 150; p <- 2; ndr <- 1
#' train <- exampleset(n, p)
#' test <- exampleset(500, p)
#'
#' # direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = 0.1, method = "direct",
#'   K = sqrt(n), keep.data = TRUE,
#'   maxitr = 150, verbose = 0, ncore = 2)
#'
#' dose <- predict(orthofit, test$X)
#' mean((test$D_opt - dose$pred)^2)
#'
#' # pseudo-direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = seq(0.1, 0.2, 0.01),
#'   method = "pseudo_direct", K = as.integer(sqrt(n)),
#'   keep.data = TRUE, maxitr = 150, ncore = 2)
#'
#' dose <- predict(orthofit, test$X)
#' mean((test$D_opt - dose$pred)^2)

orthoDr_pdose <- function(x, a, r, ndr = 2, B.initial = NULL, bw = NULL, lambda = 0.1,
                       K = sqrt(length(r)), method = c("direct","pseudo_direct"),
                       keep.data = FALSE, control = list(), maxitr = 500, verbose = 0, ncore = 0)
{
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) || nrow(x) != length(a)) stop("Number of observations do not match")

  method <- match.arg(method)
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1)
    stop("K must be a finite number greater than or equal to 1")
  K <- floor(K)
  a_scale <- sd(a)
  if (!is.finite(a_scale) || a_scale <= 0)
    stop("Observed doses must have a positive finite standard deviation")

  if (is.null(bw))
    bw = silverman(ndr, nrow(x))
  if (is.null(B.initial))
  {
    n= nrow(x)
    p = ncol(x)
    B.initial = pSAVE(x, a, r, ndr = ndr)
  }else{
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) || ndr != ncol(B.initial)) stop("Dimension of B.initial is not correct")
  }

  # check tuning parameters
  control = control.check(control)

  B.initial = gramSchmidt(B.initial)$Q

  N = nrow(x)
  P = ncol(x)
  X = x

  # standardize
  a_scale_bw = a/a_scale/bw


  cdose= seq(min(a), max(a), length.out = K)

  cdose_scale = cdose/a_scale/bw

  A.dist <- matrix(NA, nrow(x), K)
  for (k in seq_len(K))
  {
    A.dist[, k] = exp(-((a_scale_bw - cdose_scale[k]))^2)
  }

  if (method == "direct")
  {

    pre = Sys.time()
    fit = pdose_direct_solver(B.initial, X, a, A.dist, cdose, r, lambda, bw,
                           control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                           control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)
    if (verbose > 0)
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
  }

  if(method == "pseudo_direct")
  {

    pre = Sys.time()
    fit = pdose_semi_solver(B.initial, X, r, a, A.dist, cdose, lambda, bw,
                         control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                         control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)
    if (verbose > 0)
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))

  }

  fit$method = method
  fit$keep.data = keep.data

  if (keep.data)
  {
    fit[['x']] = x
    fit[['a']] = a
    fit[['r']] = r
    fit[['bw']] = bw
  }

  class(fit) <- c("orthoDr", "fit", "pdose")

  return(fit)
}

