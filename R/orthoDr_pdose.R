#' Personalized Dose Estimation via Dimension Reduction
#'
#' @description
#' Fits a personalized dose model using direct or pseudo-direct learning,
#' based on Zhou, Zhu & Zeng (2021). Direct learning targets a dose
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
#' @param lambda Positive finite penalty for kernel ridge regression. For the
#'   direct method this controls initialization only; a vector of candidates
#'   invokes GCV once at the initial direction matrix. For pseudo-direct
#'   learning it controls the final ridge dose rule. Default: \code{0.1}.
#' @param K Number of candidate dose grid points used during fitting. The grid
#'   spans the observed dose range. For direct learning, it is used only for
#'   initialization. Must be a finite scalar at least one;
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
#'   \code{gtol}. Direct learning uses analytic gradients and monotone line
#'   searches, so \code{gamma} and \code{epsilon} are unused. Its \code{gtol}
#'   is an absolute tolerance for both gradient norms; \code{btol} and
#'   \code{ftol} detect unsuccessful stagnation, not convergence.
#' @param maxitr Maximum number of iterations. Default: \code{500}. For direct
#'   learning, an iteration is a completed coefficient and direction sweep;
#'   zero returns the initialized model without optimization. Supply a finite
#'   nonnegative number within the integer range; fractional values are floored.
#' @param verbose Integer; if > 0, prints iteration progress. Default: 0.
#' @param ncore Integer; number of cores for parallel computation via OpenMP.
#'   0 = auto-detect. Default: \code{0}. The direct method's analytic solver
#'   runs serially; this argument is retained for calling compatibility.
#' @param lambda.reward Direct learning only: one finite nonnegative coefficient
#'   for the squared coefficient penalty in the reward objective. Default:
#'   \code{0.1}, a nominal choice to assess for your data. Zero removes this
#'   penalty. It is separate from \code{lambda} and is not selected by GCV.
#'   Omit this argument for \code{method = "pseudo_direct"}.
#'
#' @return An object of class \code{c("orthoDr", "fit", "pdose")}, which is a
#'   list containing:
#'   \describe{
#'     \item{\code{B}}{The estimated direction matrix (columns are orthonormal).}
#'     \item{\code{W}}{Coefficients of the continuous kernel dose rule.}
#'     \item{\code{fn}}{The final objective function value. For direct learning,
#'       this is negative estimated mean reward plus the coefficient penalty.}
#'     \item{\code{itr}}{Number of completed iterations.}
#'     \item{\code{converge}}{Logical; \code{TRUE} when a convergence criterion
#'       was reached. Direct learning requires both coefficient and constrained
#'       direction gradient norms to be at most \code{control$gtol}.}
#'     \item{\code{value}, \code{penalty}}{Direct learning: estimated mean reward
#'       and penalty contribution, with \code{fn = -value + penalty}.}
#'     \item{\code{lambda.init}, \code{lambda.reward}}{Direct learning: selected
#'       ridge initialization value and the supplied reward-penalty value.}
#'     \item{\code{gradient.norm}}{Direct learning: gradient norms named
#'       \code{B} and \code{W}, evaluated at the returned model. The direction
#'       norm uses the full Euclidean Stiefel tangent projection.}
#'     \item{\code{stop.reason}}{Direct learning: whether gradient tolerance,
#'       the iteration limit, a failed line search, stagnation, or a numerical
#'       failure stopped fitting. Only gradient tolerance indicates convergence.}
#'     \item{\code{objective.version}}{Direct learning: \code{"direct_joint_v1"},
#'       identifying the penalized joint objective.}
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
#' Direct learning initializes a kernel dose rule using grid-selected doses
#' and ridge regression. It then alternates updates of the coefficients
#' \eqn{w} and directions \eqn{B} to minimize
#' \deqn{Q(B,w) = -\widehat V(B,w) + \lambda_{\mathrm{reward}} w^\top w,}
#' subject to \eqn{B^\top B = I}. Here \eqn{\widehat V} is the estimated mean
#' reward, with no extra division of the penalty by sample size. The penalty
#' is a squared coefficient norm, not the RKHS norm \eqn{w^\top K_B w}.
#' The ridge parameter \code{lambda} is the amount added to the kernel matrix
#' during initialization and is not reused as the reward penalty.
#'
#' The dose rule uses a finite Gaussian kernel expansion centred at the
#' projected training observations. During a direction update, the coefficients
#' remain fixed while the centres and projected standard deviations change
#' with \eqn{B}. This joint finite-dimensional formulation differs from the
#' paper's literal fixed-function direction update. Ridge coefficients are not
#' refitted during optimization or after stopping.
#'
#' Pseudo-direct learning constructs a ridge dose rule after fitting the reward
#' surface. Both methods predict continuous doses, which need not lie on the
#' grid or within the observed dose range.
#'
#' Dose comparisons use the standard deviation of the observed training doses
#' for both observed and candidate doses. Returned doses remain in their
#' original units. A positive \code{lambda.reward} shrinks the dose rule toward
#' zero, so its choice depends on dose units and origin as well as reward scale.
#' Multiplying rewards by \eqn{c > 0} requires multiplying this penalty by
#' \eqn{c} to preserve the objective's tradeoff. Multiplying dose units by
#' \eqn{s > 0} requires dividing the penalty by \eqn{s^2}. There is no general
#' invariance to shifting the dose origin. The nominal default is not a
#' universally optimal choice. Assess several settings using policy evaluation
#' appropriate to your data and study design. Do not select this parameter by
#' comparing \code{fn} across settings: each uses a different penalized objective.
#'
#' Inspect \code{converge}, \code{gradient.norm}, and \code{stop.reason} before
#' using a direct fit. Failed searches preserve the last completed sweep;
#' small changes alone do not certify convergence. Convergence describes
#' numerical stationarity, not global optimality or performance on new data.
#' The returned \code{value} is an in-sample reward estimate.
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
#'   ndr = ndr, lambda = 0.1, lambda.reward = 0.1, method = "direct",
#'   K = sqrt(n), keep.data = TRUE,
#'   maxitr = 500, verbose = 0, ncore = 2)
#' orthofit$converge
#' orthofit$gradient.norm
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
                       keep.data = FALSE, control = list(), maxitr = 500, verbose = 0, ncore = 0,
                       lambda.reward = 0.1)
{
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) || nrow(x) != length(a)) stop("Number of observations do not match")

  method <- match.arg(method)
  if (method == "pseudo_direct" && !missing(lambda.reward))
    stop("lambda.reward is only used with method = 'direct'")
  if (method == "direct") {
    finite_scalar <- function(value) {
      is.numeric(value) && !is.complex(value) && length(value) == 1L &&
        is.finite(value)
    }
    if (nrow(x) < 2L || ncol(x) < 1L || is.complex(x) || any(!is.finite(x)))
      stop("x must have at least two rows and one column, with finite numeric values")
    for (name in c("a", "r")) {
      value <- get(name)
      if (!is.numeric(value) || is.complex(value) || !is.null(dim(value)) ||
          any(!is.finite(value)))
        stop(name, " must be a finite numeric vector")
    }
    if (!finite_scalar(ndr) || ndr < 1 || ndr != floor(ndr) || ndr > ncol(x))
      stop("ndr must be an integer between 1 and ncol(x)")
    if (!is.numeric(lambda) || is.complex(lambda) || length(lambda) == 0L ||
        !is.null(dim(lambda)) || any(!is.finite(lambda)) || any(lambda <= 0))
      stop("lambda must contain positive finite ridge parameters")
    if (!finite_scalar(lambda.reward) || lambda.reward < 0)
      stop("lambda.reward must be one finite nonnegative number")
    if (!finite_scalar(maxitr) || maxitr < 0 || floor(maxitr) > .Machine$integer.max)
      stop("maxitr must be a finite nonnegative number within the integer range")
    maxitr <- as.integer(floor(maxitr))
    if (!finite_scalar(ncore) || ncore < 0 || ncore != floor(ncore) ||
        ncore > .Machine$integer.max)
      stop("ncore must be a nonnegative integer within the integer range")
    if (is.logical(verbose) && length(verbose) == 1L && !is.na(verbose))
      verbose <- as.integer(verbose)
    if (!finite_scalar(verbose) || verbose < 0 || verbose != floor(verbose) ||
        verbose > .Machine$integer.max)
      stop("verbose must be a nonnegative integer within the integer range")
    if (!is.list(control)) stop("control must be a list of tuning parameters")
    defaults <- list(rho = 1e-4, eta = 0.2, tau = 1e-3,
                     btol = 1e-6, ftol = 1e-6, gtol = 1e-6)
    for (name in names(defaults)) {
      if (!name %in% names(control)) control[[name]] <- defaults[[name]]
      value <- control[[name]]
      if (!finite_scalar(value) || value <= 0 ||
          (name %in% c("rho", "eta") && value >= 1))
        stop("control$", name, " must be a positive finite scalar",
             if (name %in% c("rho", "eta")) " less than 1" else "")
    }
    # These arguments remain in the native interface for other callers.
    # Direct learning now uses analytic gradients and monotone line searches.
    control$gamma <- 0.85
    control$epsilon <- 1e-6
  }
  if (method == "direct" &&
      (!finite_scalar(K) || K < 1 || floor(K) > .Machine$integer.max))
    stop("K must be a finite real number at least 1 within the integer range")
  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) || K < 1)
    stop("K must be a finite number greater than or equal to 1")
  K <- floor(K)
  a_scale <- sd(a)
  if (!is.finite(a_scale) || a_scale <= 0)
    stop("Observed doses must have a positive finite standard deviation")

  if (is.null(bw))
    bw = silverman(ndr, nrow(x))
  if (method == "direct" && (!finite_scalar(bw) || bw <= 0))
    stop("bw must be one positive finite number")
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
  if (method == "pseudo_direct") control = control.check(control)

  if (method == "direct" &&
      (!is.numeric(B.initial) || is.complex(B.initial) || any(!is.finite(B.initial)) ||
       qr(B.initial)$rank < ndr))
    stop("B.initial must be finite numeric with linearly independent columns")

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
                           control$btol, control$ftol, control$gtol, maxitr, verbose, ncore,
                           lambda.reward)
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

