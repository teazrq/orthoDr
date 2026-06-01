#' Orthogonality Constrained Optimization
#'
#' @description
#' A general-purpose optimization solver for problems with orthogonality
#' constraints on the columns of \code{B}. Implements the curvilinear search
#' algorithm of Wen & Yin (2013), which ensures that \code{B} remains on the
#' Stiefel manifold \eqn{B^\top B = I} throughout optimization.
#' If no analytical gradient is provided, a numerical approximation is used.
#'
#' @name ortho_optim
#' @export
#'
#' @param B Initial \code{B} values. A matrix whose columns are subject to the
#'   orthogonality constraint. If the columns are not orthogonal, they will be
#'   orthogonalized via Gram-Schmidt before optimization begins.
#' @param fn Objective function. The first argument must be \code{B} (a matrix).
#'   Should return a single numeric value.
#' @param grad Gradient function. The first argument must be \code{B}. Should
#'   return a matrix of the same dimension as \code{B}. If \code{NULL} (default),
#'   the gradient is approximated numerically via finite differences.
#' @param ... Additional named arguments passed to \code{fn} and \code{grad}.
#' @param maximize Logical; if \code{TRUE}, the solver maximizes \code{fn}
#'   instead of minimizing it. Default is \code{FALSE}.
#' @param control A list of tuning parameters for the optimizer:
#'   \describe{
#'     \item{\code{rho}}{Controls the linear approximation quality in the
#'       curvilinear search. Default: \code{1e-4}.}
#'     \item{\code{eta}}{Factor for decreasing the step size in backtracking
#'       line search. Default: \code{0.2}.}
#'     \item{\code{gamma}}{Nonmonotone search parameter from Zhang & Hager
#'       (2004). Default: \code{0.85}.}
#'     \item{\code{tau}}{Initial step size for updating \code{B}.
#'       Default: \code{1e-3}.}
#'     \item{\code{epsilon}}{Step size for finite-difference gradient
#'       approximation (only used when \code{grad = NULL}).
#'       Default: \code{1e-6}.}
#'     \item{\code{btol}}{Tolerance for orthogonality deviation of \code{B}.
#'       Convergence is reached when \eqn{\|B^\top B - I\|_F < \code{btol}}.
#'       Default: \code{1e-6}.}
#'     \item{\code{ftol}}{Tolerance for relative change in the objective value.
#'       Default: \code{1e-6}.}
#'     \item{\code{gtol}}{Tolerance for the gradient (projected onto the
#'       manifold). Default: \code{1e-6}.}
#'   }
#' @param maxitr Maximum number of iterations. Default: \code{500}.
#' @param verbose Integer; if > 0, prints iteration progress. Default: 0.
#'
#' @return An object of class \code{c("orthoDr", "fit", "optim")}, which is a
#'   list containing:
#'   \describe{
#'     \item{\code{B}}{The optimized matrix (columns are orthonormal).}
#'     \item{\code{fn}}{The final objective function value.}
#'     \item{\code{fn_Seq}}{Numeric vector of objective values at each iteration
#'       (length \code{maxitr}, padded with zeros after convergence).}
#'     \item{\code{itr}}{Number of iterations performed.}
#'     \item{\code{converge}}{Convergence code.}
#'     \item{\code{method}}{\code{"true gradient"} or \code{"approx. gradient"}.}
#'   }
#'
#' @references
#' Wen, Z. and Yin, W. (2013). A feasible method for optimization with
#' orthogonality constraints. \emph{Mathematical Programming}, 142(1-2),
#' 397--434. \doi{10.1007/s10107-012-0584-1}.
#'
#' Zhang, H. and Hager, W. W. (2004). A nonmonotone line search technique
#' and its application to unconstrained optimization.
#' \emph{SIAM Journal on Optimization}, 14(4), 1043--1056.
#' \doi{10.1137/S1052623403426556}.
#'
#' @examples
#' # Eigenvalue problem: minimize -0.5 * tr(B'A B) s.t. B'B = I
#' library(pracma)
#' set.seed(1)
#' n <- 100; k <- 6
#' A <- matrix(rnorm(n * n), n, n)
#' A <- t(A) %*% A
#' B <- gramSchmidt(matrix(rnorm(n * k), n, k))$Q
#'
#' fx <- function(B, A) -0.5 * sum(diag(t(B) %*% A %*% B))
#' gx <- function(B, A) -A %*% B
#' fit <- ortho_optim(B, fx, gx, A = A)
#' fx(fit$B, A)
#'
#' # Compare with the analytical solution from eigen()
#' sol <- eigen(A)$vectors[, 1:k]
#' fx(sol, A)


ortho_optim <- function(B, fn, grad = NULL, ..., maximize = FALSE,
                        control = list(), maxitr = 500, verbose = 0)
{
  if (is.null(B))
  {
    stop("Initial value of B must be given")
  }else{
    if (any(is.na(B))) stop("B cannot contain NA values")
    if (!is.matrix(B)) stop("B must be a matrix")
  }

  # check orthogonality of initial value
  if (sum(abs(t(B) %*% B - diag(ncol(B)))) > 1e-15)
  {
    cat("Initial B not orthogonal, will be processed by Gram-Schmidt \n")
    B = gramSchmidt(B)$Q
  }

  if (verbose)
  {
    cat(paste("Optimizing", dim(B)[1]*dim(B)[2], "parameters with", dim(B)[2], "orthogonality constraints...\n"))
  }

  # check tuning parameters

  control = control.check(control)

  # check objects

  env = environment()

  names = sapply(substitute(list(...))[-1], deparse)

  if (length(names) > 0)
  {
    for (i in 1:length(names))
      if (!exists(names[[i]], envir = env))
        stop(paste(names[[i]], "do not exist"))
  }

  # check f and g

  if (is.null(fn))
    stop("fn must be given")

  if (maximize)
    f <- function(par) -fn(par, ...)
  else
    f <- function(par) fn(par, ...)

  fB <- f(B)
  if (!is.numeric(fB) || length(fB) != 1)
    stop("fn must return a single number")

  if (!is.null(grad))
  {
    useg = TRUE

    if (maximize)
      g <- function(par) -grad(par, ...)
    else
      g <- function(par) grad(par, ...)

    gB <- g(B)
    if (!is.matrix(gB))
      stop("grad must return a matrix")
    else if (any(dim(gB) != dim(B)))
      stop("grad must return a matrix with the same dimension as B")
  }else{
    useg = FALSE
    g <- function(par) stop("cannot use grad")
  }

  pre = Sys.time()
  fit = gen_solver(B, f, g, env, useg, control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                   control$btol, control$ftol, control$gtol, maxitr, verbose)
  if (verbose > 0)
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))

  if (maximize)
    fit$fn = -fit$fn

  fit$method = ifelse(useg, "true gradient", "approx. gradient")

  class(fit) <- c("orthoDr", "fit", "optim")

  return(fit)
}
