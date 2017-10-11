#' @title Orthogonality constrained optimization
#' @name ortho_optim
#' @description A general purpose optimization solver with orthogonality constraint. The orthogonality constrained optimization method is a nearly direct translation from Wen and Yin (2010)'s matlab code.
#' @param B Initial \code{B} values. Must be a matrix, and the columns are subject to the orthogonality constrains. Will be processed by Gram-Schmidt if not orthogonal
#' @param fn A function that calculate the objective function value. The first argument should be \code{B}. Returns a single value.
#' @param grad A function that calculate the gradient. The first argument should be \code{B}. Returns a matrix with the same dimension as \code{B}. If not specified, then numerical approximation is used.
#' @param ... Arguments passed to \code{fn} and \code{grad}
#' @param maximize By default, the solver will try to minimize the objective function unless \code{maximize = TRUE}
#' @param control A list of tuning varaibles for optimization. \code{epsilon} is the size for numerically appriximating the gradient. For others, see Wen and Yin (2013).
#' @param maxitr Maximum number of iterations
#' @param verbose Should information be displayed
#' @return A \code{orthoDr} object; a list consisting of
#' \item{B}{The optimal \code{B} value}
#' \item{fn}{The final funtional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434.
#' DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
#' @examples
#' library(pracma)
#' # an example of searching for the first principal component
#' set.seed(1)
#' P = 100
#' N = 400
#' X = scale(matrix(rnorm(N*P), N, P), scale = FALSE)
#' w = gramSchmidt(matrix(rnorm(P), P, 1))$Q
#' fx <- function(w, X) t(w) %*% t(X) %*% X %*% w
#'
#' # fit the
#' fit = ortho_optim(w, fx, X = X, maximize = TRUE, verbose = 0)
#'
#' # compare the result with prcomp
#' distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist")
#'
#' # if we have the derivative
#'
#' gx <- function(w, X) 2*t(X) %*% X %*% w
#' fit = ortho_optim(w, fx, gx, X = X, maximize = TRUE, verbose = 0)
#' distance(fit$B, as.matrix(prcomp(X)$rotation[, 1]), type = "dist")

ortho_optim <- function(B, fn, grad = NULL, ..., maximize = FALSE,
                        control = list(), maxitr = 500, verbose = FALSE)
{
  if (is.null(B))
  {
    stop("Initial value of B must be given")
  }else{
    if (any(is.na(B))) stop("B cannot contain NA values")
    if (!is.matrix(B)) stop("B must be a matrix")
  }

  # center matrix X, but do not scale
  if (sum(abs(t(B) %*% B - diag(ncol(B)))) > 1e-15)
  {
    cat("Initial B not orthognal, will be processed by Gram-Schmidt \n")
    B = gramSchmidt(B)$Q
  }

  if (verbose)
  {
    cat(paste("Optimizing", dim(B)[1]*dim(B)[2], "parameters with", dim(B)[2], "orthogonality constrains...\n"))
  }

  # check tuning parameters

  if (!is.list(control))
    stop("control must be a list of tuning parameters")

  if (is.null(control$rho))
  {
    control$rho = 1e-4
  }else if(control$rho < 0) control$rho = 1e-4

  if (is.null(control$eta))
  {
    control$eta = 0.2
  }else if(control$eta < 0) control$eta = 0.2

  if (is.null(control$gamma))
  {
    control$gamma = 0.85
  }else if(control$gamma < 0) control$gamma = 0.85

  if (is.null(control$tau))
  {
    control$tau = 1e-3
  }else if(control$tau < 0) control$tau = 1e-3

  if (is.null(control$epsilon))
  {
    control$epsilon = 1e-6
  }else if(control$epsilon < 0) control$epsilon = 1e-6

  if (is.null(control$btol))
  {
    control$btol = 1e-6
  }else if(control$btol < 0) control$btol = 1e-6

  if (is.null(control$ftol))
  {
    control$ftol = 1e-6
  }else if(control$ftol < 0) control$ftol = 1e-6

  if (is.null(control$gtol))
  {
    control$gtol = 1e-6
  }else if(control$gtol < 0) control$gtol = 1e-6

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

  if ( !is.numeric(f(B)) |  (length(f(B)) != 1) )
    stop("fn must return a single number")

  if (!is.null(grad))
  {
    useg = TRUE

    if (maximize)
      g <- function(par) -grad(par, ...)
    else
      g <- function(par) grad(par, ...)

    if (!is.matrix(g(B)))
      stop("grad must return a matrix")
    else if(  any(dim(g(B)) != dim(B)) )
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

  class(fit) <- "orthoDr"

  return(fit)
}
