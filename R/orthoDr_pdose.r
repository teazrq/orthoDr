#' Direct Learning & Pseudo-direct Learning Model
#' 
#' Performs the "Direct Learning & Pseudo-direct Learning" Method for
#' personalized medicine.
#' 
#' @param x         A `matrix` or `data.frame` for features (continuous only).
#' @param a         A `vector` of observed dose
#' @param r         A `vector` of observed reward
#' @param ndr       A dimension structure
#' @param B.initial Initial `B` values. Will use the partial SAVE [pSAVE][orthoDr::pSAVE] as the initial if 
#'                  leaving as `NULL`. If specified, must be a matrix with 
#'                  `ncol(x)` rows and `ndr` columns. Will be processed by
#'                   Gram-Schmidt if not orthogonal.
#' @param bw        A Kernel bandwidth, assuming each variables have unit variance
#' @param lambda    The penalty level for kernel ridge regression. If a range of values is specified, the GCV will be used to select the best tuning
#' @param K         A number of grids in the range of dose
#' @param method    Either `"direct"` or `"pseudo_direct"` 
#' @param keep.data Should the original data be kept for prediction
#' @param control   A list of tuning variables for optimization. `epsilon` is the size for numerically approximating the gradient. For others, see Wen and Yin (2013).
#' @param maxitr    Maximum number of iterations
#' @param ncore     the number of cores for parallel computing
#' @param verbose   Should information be displayed
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
#' Zhou, W., Zhu, R., & Zeng, D. (2021). A parsimonious personalized dose-finding model via dimension reduction. 
#' Biometrika, 108(3), 643-659.
#' DOI: \doi{https://doi.org/10.1093/biomet/asaa087}
#' 
#' @examples
#' # generate some personalized dose scenario
#'
#' exampleset <- function(size, ncov) {
#'   X <- matrix(runif(size * ncov, -1, 1), ncol = ncov)
#'   A <- runif(size, 0, 2)
#'
#'   Edr <- as.matrix(c(0.5, -0.5))
#'
#'   D_opt <- X %*% Edr + 1
#'
#'   mu <- 2 + 0.5 * (X %*% Edr) - 7 * abs(D_opt - A)
#'
#'   R <- rnorm(length(mu), mu, 1)
#'
#'   R <- R - min(R)
#'
#'   datainfo <- list(X = X, A = A, R = R, D_opt = D_opt, mu = mu)
#'   return(datainfo)
#' }
#'
#' # generate data
#'
#' set.seed(123)
#' n <- 150
#' p <- 2
#' ndr <- 1
#' train <- exampleset(n, p)
#' test <- exampleset(500, p)
#'
#' # the direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = 0.1,
#'   method = "direct", K = sqrt(n), keep.data = TRUE,
#'   maxitr = 150, verbose = FALSE, ncore = 2
#' )
#'
#' dose <- predict(orthofit, test$X)
#'
#' # ` # compare with the optimal dose
#' dosedistance <- mean((test$D_opt - dose$pred)^2)
#' print(dosedistance)
#'
#' # the pseudo direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = seq(0.1, 0.2, 0.01),
#'   method = "pseudo_direct", K = as.integer(sqrt(n)), keep.data = TRUE,
#'   maxitr = 150, verbose = FALSE, ncore = 2
#' )
#'
#' dose <- predict(orthofit, test$X)
#'
#' # compare with the optimal dose
#'
#' dosedistance <- mean((test$D_opt - dose$pred)^2)
#' print(dosedistance)
orthoDr_pdose <- function(x, a, r, ndr = ndr, B.initial = NULL, bw = NULL, lambda = 0.1,
                          K = sqrt(length(r)), method = c("direct", "pseudo_direct"),
                          keep.data = FALSE, control = list(), maxitr = 500, verbose = FALSE, ncore = 0) {
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) | nrow(x) != length(a)) stop("Number of observations do not match")

  if (is.null(bw)) {
    bw <- silverman(ndr, nrow(x))
  }
  if (is.null(B.initial)) {
    n <- nrow(x)
    p <- ncol(x)
    B.initial <- pSAVE(x, a, r, ndr = ndr)
  } else {
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) | ndr != ncol(B.initial)) stop("Dimention of B.initial is not correct")
  }

  # check tuning parameters
  control <- control.check(control)

  B.initial <- gramSchmidt(B.initial)$Q

  N <- nrow(x)
  P <- ncol(x)
  X <- x

  # standerdize
  a_center <- mean(a)
  a_scale <- sd(a)
  a_scale_bw <- a / a_scale / bw


  cdose <- seq(min(a), max(a), length.out = K)

  cdose_scale <- cdose / sd(cdose) / bw

  A.dist <- matrix(NA, nrow(x), K)
  for (k in 1:K)
  {
    A.dist[, k] <- exp(-((a_scale_bw - cdose_scale[k]))^2)
  }

  if (method == "direct") {
    pre <- Sys.time()
    fit <- pdose_direct_solver(
      B.initial, X, a, A.dist, cdose, r, lambda, bw,
      control$rho, control$eta, control$gamma, control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
    )
    if (verbose > 0) {
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
    }
  }

  if (method == "pseudo_direct") {
    pre <- Sys.time()
    fit <- pdose_semi_solver(
      B.initial, X, r, a, A.dist, cdose, lambda, bw,
      control$rho, control$eta, control$gamma, control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
    )
    if (verbose > 0) {
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
    }
  }

  fit$method <- method
  fit$keep.data <- keep.data

  if (keep.data) {
    fit[["x"]] <- x
    fit[["a"]] <- a
    fit[["r"]] <- r
    fit[["bw"]] <- bw
  }

  class(fit) <- c("orthoDr", "fit", "pdose")

  return(fit)
}
