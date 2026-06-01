#' Semiparametric Dimension Reduction for Censored Survival Outcomes
#'
#' @description
#' Fits a counting-process-based semiparametric dimension reduction (IR-CP)
#' model for right-censored survival outcomes (Sun, Zhu, Wang & Zeng, 2017).
#' Three estimating equations are available: \code{"forward"}, \code{"dn"},
#' and \code{"dm"}.
#'
#' @name orthoDr_surv
#' @export
#'
#' @param x A matrix or data.frame of features. Columns are not automatically
#'   scaled to unit variance.
#' @param y A numeric vector of observed survival times.
#' @param censor A numeric vector of censoring indicators (1 = event,
#'   0 = censored).
#' @param method The estimating equation: \code{"forward"} (1-d model),
#'   \code{"dn"} (counting process), or \code{"dm"} (martingale).
#'   Default: \code{"dm"}.
#' @param ndr Number of directions to estimate. Default: \code{2}
#'   (automatically set to \code{1} for \code{method = "forward"}).
#' @param B.initial Initial \code{B} matrix. If \code{NULL} (default),
#'   the counting-process-based SIR model (\code{\link{CP_SIR}}) is used.
#'   If specified, must have \code{nrow(B.initial) == ncol(x)} and
#'   \code{ncol(B.initial) == ndr}. Will be orthogonalized via Gram-Schmidt.
#' @param bw Kernel bandwidth, assuming each variable has unit variance.
#'   If \code{NULL} (default), Silverman's rule is used.
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
#' @return An object of class \code{c("orthoDr", "fit", "surv")}, which is a
#'   list containing:
#'   \describe{
#'     \item{\code{B}}{The estimated direction matrix (columns are orthonormal).}
#'     \item{\code{fn}}{The final objective function value.}
#'     \item{\code{itr}}{Number of iterations performed.}
#'     \item{\code{converge}}{Convergence code.}
#'     \item{\code{method}}{The method used (e.g., \code{"dm"}).}
#'     \item{\code{keep.data}}{Whether original data was retained.}
#'   }
#'
#' @references
#' Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2017). Counting process based
#' dimension reduction method for censored outcomes.
#' \emph{arXiv preprint} \doi{10.48550/arXiv.1704.05046}.
#'
#' @examples
#' # Setting 1 from Sun et al. (2017), reduced sample size
#' library(MASS)
#' set.seed(1)
#' N <- 200; P <- 6
#' V <- 0.5^abs(outer(1:P, 1:P, "-"))
#' dataX <- as.matrix(mvrnorm(N, mu = rep(0, P), Sigma = V))
#' failEDR <- as.matrix(c(1, 0.5, 0, 0, 0, rep(0, P - 5)))
#' censorEDR <- as.matrix(c(0, 0, 0, 1, 1, rep(0, P - 5)))
#' T <- rexp(N, exp(dataX %*% failEDR))
#' C <- rexp(N, exp(dataX %*% censorEDR - 1))
#' ndr <- 1
#' Y <- pmin(T, C)
#' Censor <- (T < C)
#'
#' # forward model (1-d only)
#' forward.fit <- orthoDr_surv(dataX, Y, Censor, method = "forward")
#' distance(failEDR, forward.fit$B, "dist")
#'
#' # counting process model
#' dn.fit <- orthoDr_surv(dataX, Y, Censor, method = "dn", ndr = ndr)
#' distance(failEDR, dn.fit$B, "dist")
#'
#' # martingale model
#' dm.fit <- orthoDr_surv(dataX, Y, Censor, method = "dm", ndr = ndr)
#' distance(failEDR, dm.fit$B, "dist")
#'
#' @seealso \code{\link{orthoDr_pdose}}, \code{\link{CP_SIR}}, \code{\link{view_dr_surv}}

orthoDr_surv <- function(x, y, censor, method = "dm", ndr = ifelse(method == "forward", 1, 2),
                         B.initial = NULL,
                         bw = NULL,
                         keep.data = FALSE,
                         control = list(),
                         maxitr = 500,
                         verbose = 0,
                         ncore = 0)
{
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y) || nrow(x) != length(censor)) stop("Number of observations do not match")

  # check tuning parameters
  control = control.check(control)
  method = match.arg(method, c("forward", "dn", "dm"))

  ndr = max(1, ndr)
  if (ndr > 4) warning("ndr > 4 is not recommended due to nonparametric kernel estimations")
  ndr = min(ndr, ncol(x))

  if (method == "forward" && ndr > 1)
  {
    warning("forward can only solve for 1 dimension")
    ndr = 1
  }

  N = nrow(x)
  P = ncol(x)

  if (is.null(bw))
    bw = silverman(ndr, N)

  # scale X and Y for stability
  # Y = scale(log(y)) / sqrt(2) / silverman(1, N)

  Y = scale(rank(y, ties.method = "average") / (N + 1)) / sqrt(2) / silverman(1, N)

  xscale = apply(x, 2, sd)
  X = scale(x)

  # get initial value
  if (is.null(B.initial))
  {
    B.initial = CP_SIR(X, Y, censor)$vectors[,1:ndr, drop = FALSE]
  }else{
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) || ndr != ncol(B.initial)) stop("Dimension of B.initial is not correct")
    if (method == "forward" && ncol(B.initial) > 1) stop("forward method can only use 1-d B.initial")
  }

  B.initial = gramSchmidt(B.initial)$Q

  # pre-process

  Yorder = order(Y)
  X = X[Yorder, ]
  Y = Y[Yorder]
  C = censor[Yorder]
  Fail.Ind = which(C==1)

  # calculate some useful stuff

  nFail = length(Fail.Ind)

  # E[X | dN(t) = 1, Y(t) = 1] - E[X | dN(t) = 0, Y(t) = 1] at all failure times

  kernel.y = exp(-(as.matrix(dist(Y, method = "euclidean"))))
  kernel.y[upper.tri(kernel.y)] = 0

  Phit = matrix(0, P, nFail)

  for (j in 1:nFail)
    Phit[,j] = as.matrix(apply(X, 2, weighted.mean, w = C*kernel.y[, Fail.Ind[j]]) - colMeans(X[Fail.Ind[j]:N, ,drop = FALSE]))

  # start to fit the model

  pre = Sys.time()

  if (method == "forward")
    fit = surv_forward_solver(B.initial, X, Fail.Ind, bw,
                         control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                         control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)

  if (method == "dn")
    fit = surv_dn_solver(B.initial, X, Phit, Fail.Ind, bw,
                    control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                    control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)

  if (method == "dm")
    fit = surv_dm_solver(B.initial, X, Phit, Fail.Ind, bw,
                         control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                         control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)

  if (verbose > 0)
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))

  # rescale B back to the original scale

  fit$B = sweep(fit$B, 1, xscale, FUN = "/")
  fit$B = apply(fit$B, 2, function(x) x / sqrt(sum(x^2)))
  fit$method = method
  fit$keep.data = keep.data

  if (keep.data)
  {
    fit[['x']] = x
    fit[['y']] = y
    fit[['censor']] = censor
  }

  class(fit) <- c("orthoDr", "fit", "surv")

  return(fit)
}
