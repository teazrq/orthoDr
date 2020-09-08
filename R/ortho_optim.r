#' Orthogonality constrained optimization
#' 
#' A general purpose optimization solver with orthogonality constraint.
#' The orthogonality constrained optimization method is a nearly direct
#' translation from Wen and Yin (2010)'s `MATLAB` code.
#' 
#' @param B        Initial `B` values. Must be a matrix, and the columns are 
#'                 subject to the orthogonality constrains. Will be processed
#'                 by Gram-Schmidt if not orthogonal
#' @param fn       A function that calculate the objective function value. 
#'                 The first argument should be `B`. Returns a single value.
#' @param grad     A function that calculate the gradient. The first argument
#'                 should be `B`. Returns a matrix with the same dimension
#'                 as `B`. If not specified, then numerical approximation is
#'                 used.
#' @param ...      Arguments passed to `fn` and `grad`
#' @param maximize By default, the solver will try to minimize the objective
#'                 function unless `maximize = TRUE`
#' @param control  A list of tuning variables for optimization. `epsilon` is
#'                 the size for numerically approximating the gradient.
#'                 For others, see Wen and Yin (2013).
#' @param maxitr   Maximum number of iterations
#' @param verbose  Should information be displayed
#' 
#' @return 
#' 
#' A `orthoDr` object that consists of a `list` with named entries of:
#' 
#' \item{B}{The optimal `B` value}
#' \item{fn}{The final functional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' 
#' @export
#' 
#' @references
#' Wen, Z., & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. Mathematical Programming, 142(1), 397-434.
#' DOI: \doi{10.1007/s10107-012-0584-1}
#'
#' @examples
#' # an eigen value problem
#' library(pracma)
#' set.seed(1)
#' n <- 100
#' k <- 6
#' A <- matrix(rnorm(n * n), n, n)
#' A <- t(A) %*% A
#' B <- gramSchmidt(matrix(rnorm(n * k), n, k))$Q
#'
#' fx <- function(B, A) -0.5 * sum(diag(t(B) %*% A %*% B))
#' gx <- function(B, A) -A %*% B
#' fit <- ortho_optim(B, fx, gx, A = A)
#' fx(fit$B, A)
#'
#' # compare with the solution from the eigen function
#' sol <- eigen(A)$vectors[, 1:k]
#' fx(sol, A)
ortho_optim <- function(B, fn, grad = NULL, ..., maximize = FALSE,
                        control = list(), maxitr = 500, verbose = FALSE) {
  if (is.null(B)) {
    stop("Initial value of B must be given")
  } else {
    if (any(is.na(B))) stop("B cannot contain NA values")
    if (!is.matrix(B)) stop("B must be a matrix")
  }

  # check orthogonality of initial value
  if (sum(abs(t(B) %*% B - diag(ncol(B)))) > 1e-15) {
    cat("Initial B not orthogonal, will be processed by Gram-Schmidt \n")
    B <- gramSchmidt(B)$Q
  }

  if (verbose) {
    cat(paste("Optimizing", dim(B)[1] * dim(B)[2], "parameters with", dim(B)[2], "orthogonality constrains...\n"))
  }

  # check tuning parameters

  control <- control.check(control)

  env <- environment()

  # check f and g

  if (is.null(fn)) {
    stop("fn must be given")
  }

  if (maximize) {
    f <- function(par) -fn(par, ...)
  } else {
    f <- function(par) fn(par, ...)
  }

  if (!is.numeric(f(B)) | (length(f(B)) != 1)) {
    stop("fn must return a single number")
  }

  if (!is.null(grad)) {
    useg <- TRUE

    if (maximize) {
      g <- function(par) -grad(par, ...)
    } else {
      g <- function(par) grad(par, ...)
    }

    if (!is.matrix(g(B))) {
      stop("grad must return a matrix")
    } else if (any(dim(g(B)) != dim(B))) {
      stop("grad must return a matrix with the same dimension as B")
    }
  } else {
    useg <- FALSE
    g <- function(par) stop("cannot use grad")
  }

  pre <- Sys.time()
  fit <- gen_solver(
    B, f, g, env, useg, control$rho, control$eta, control$gamma, control$tau, control$epsilon,
    control$btol, control$ftol, control$gtol, maxitr, verbose
  )
  if (verbose > 0) {
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
  }

  if (maximize) {
    fit$fn <- -fit$fn
  }

  fit$method <- ifelse(useg, "true gradient", "approx. gradient")

  class(fit) <- c("orthoDr", "fit", "optim")

  return(fit)
}
