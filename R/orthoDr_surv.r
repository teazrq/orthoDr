#' Counting Process based semiparametric dimension reduction (IR-CP) model
#' 
#' Models the data according to the counting process based semiparametric
#' dimension reduction (IR-CP) model for right censored survival outcome.
#' 
#' @param x         A `matrix` or `data.frame` for features. 
#'                  The algorithm will not scale the columns to unit variance
#' @param y         A `vector` of observed time
#' @param censor    A `vector` of censoring indicator
#' @param method    Estimation equation to use. Either:
#'                  `"forward"` (1-d model), `"dn"` (counting process), or
#'                   `"dm"` (martingale).
#' @param ndr       The number of directions
#' @param B.initial Initial `B` values. Will use the counting process based
#'                  SIR model [CP_SIR][orthoDr::CP_SIR] as the initial if
#'                  leaving as `NULL`. If specified, must be a matrix with 
#'                  `ncol(x)` rows and `ndr` columns. Will be processed by
#'                  Gram-Schmidt if not orthogonal.
#' @param bw        A Kernel bandwidth, assuming each variables have
#'                  unit variance.
#' @param keep.data Should the original data be kept for prediction. 
#'                  Default is `FALSE`
#' @param control   A list of tuning variables for optimization. `epsilon` is
#'                  the size for numerically approximating the gradient. 
#'                  For others, see Wen and Yin (2013).
#' @param maxitr    Maximum number of iterations
#' @param verbose   Should information be displayed
#' @param ncore     Number of cores for parallel computing.
#'                  The default is the maximum number of threads.
#' 
#' @return 
#' 
#' A `orthoDr` object consisting of `list` with named elements:
#' 
#' \item{B}{The optimal `B` value}
#' \item{fn}{The final functional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' 
#' @export
#' 
#' @references 
#' Sun, Q., Zhu, R., Wang, T., & Zeng, D. (2019). Counting process-based dimension reduction methods for censored outcomes. Biometrika, 106(1), 181-196.
#' DOI: \doi{10.1093/biomet/asy064}
#' 
#' @examples
#' # This is setting 1 in Sun et. al. (2017) with reduced sample size
#' library(MASS)
#' set.seed(1)
#' N <- 200
#' P <- 6
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
#' # fit the model
#' forward.fit <- orthoDr_surv(dataX, Y, Censor, method = "forward")
#' distance(failEDR, forward.fit$B, "dist")
#'
#' dn.fit <- orthoDr_surv(dataX, Y, Censor, method = "dn", ndr = ndr)
#' distance(failEDR, dn.fit$B, "dist")
#'
#' dm.fit <- orthoDr_surv(dataX, Y, Censor, method = "dm", ndr = ndr)
#' distance(failEDR, dm.fit$B, "dist")
orthoDr_surv <- function(x, y, censor, method = "dm", ndr = ifelse(method == "forward", 1, 2),
                         B.initial = NULL,
                         bw = NULL,
                         keep.data = FALSE,
                         control = list(),
                         maxitr = 500,
                         verbose = FALSE,
                         ncore = 0) {
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y) | nrow(x) != length(censor)) stop("Number of observations do not match")

  # check tuning parameters
  control <- control.check(control)
  match.arg(method, c("forward", "dn", "dm"))

  ndr <- max(1, ndr)
  if (ndr > 4) warning("ndr > 4 is not recommended due to nonparametric kernel estimations")
  ndr <- min(ndr, ncol(x))

  if (method == "forward" & ndr > 1) {
    warning("forward can only solve for 1 dimension")
    ndr <- 1
  }

  N <- nrow(x)
  P <- ncol(x)

  if (is.null(bw)) {
    bw <- silverman(ndr, N)
  }

  # scale X and Y for stability
  # Y = scale(log(y)) / sqrt(2) / silverman(1, N)

  Y <- scale(rank(y, ties.method = "average") / (N + 1)) / sqrt(2) / silverman(1, N)

  xscale <- apply(x, 2, sd)
  X <- scale(x)

  # get initial value
  if (is.null(B.initial)) {
    B.initial <- CP_SIR(X, Y, censor)$vectors[, 1:ndr, drop = FALSE]
  } else {
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) | ndr != ncol(B.initial)) stop("Dimension of B.initial is not correct")
    if (method == "forward" & ncol(B.initial) > 1) stop("forward method can only use 1-d B.initial")
  }

  B.initial <- gramSchmidt(B.initial)$Q

  # pre-process

  Yorder <- order(Y)
  X <- X[Yorder, ]
  Y <- Y[Yorder]
  C <- censor[Yorder]
  Fail.Ind <- which(C == 1)

  # calculate some useful stuff

  nFail <- length(Fail.Ind)

  # E[X | dN(t) = 1, Y(t) = 1] - E[X | dN(t) = 0, Y(t) = 1] at all failure times

  kernel.y <- exp(-(as.matrix(dist(Y, method = "euclidean"))))
  kernel.y[upper.tri(kernel.y)] <- 0

  Phit <- matrix(0, P, nFail)

  for (j in 1:nFail) {
    Phit[, j] <- as.matrix(apply(X, 2, weighted.mean, w = C * kernel.y[, Fail.Ind[j]]) - colMeans(X[Fail.Ind[j]:N, , drop = FALSE]))
  }

  # start to fit the model

  pre <- Sys.time()

  if (method == "forward") {
    fit <- surv_forward_solver(
      B.initial, X, Fail.Ind, bw,
      control$rho, control$eta, control$gamma, control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
    )
  }

  if (method == "dn") {
    fit <- surv_dn_solver(
      B.initial, X, Phit, Fail.Ind, bw,
      control$rho, control$eta, control$gamma, control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
    )
  }

  if (method == "dm") {
    fit <- surv_dm_solver(
      B.initial, X, Phit, Fail.Ind, bw,
      control$rho, control$eta, control$gamma, control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
    )
  }

  if (verbose > 0) {
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
  }

  # rescale B back to the original scale

  fit$B <- sweep(fit$B, 1, xscale, FUN = "/")
  fit$B <- apply(fit$B, 2, function(x) x / sqrt(sum(x^2)))
  fit$method <- method
  fit$keep.data <- keep.data

  if (keep.data) {
    fit[["x"]] <- x
    fit[["y"]] <- y
    fit[["censor"]] <- censor
  }

  class(fit) <- c("orthoDr", "fit", "surv")

  return(fit)
}
