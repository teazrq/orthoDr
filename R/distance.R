#' @title Distance Between Two Linear Subspaces
#' @name distance
#' @export
#'
#' @description
#' Compares the linear subspaces spanned by the columns of \code{s1} and
#' \code{s2} using projection matrices or correlations.
#'
#' @param s1 A matrix with linearly independent columns spanning the first subspace.
#' @param s2 A matrix with linearly independent columns spanning the second subspace,
#'   with the same number of rows as \code{s1}.
#' @param type Character string specifying the comparison to use.
#'   \describe{
#'     \item{\code{"dist"}}{Frobenius norm of the difference between
#'       projection matrices (default).}
#'     \item{\code{"trace"}}{Trace correlation. For subspaces of equal
#'       dimension, 1 indicates identical subspaces.}
#'     \item{\code{"canonical"}}{Mean canonical correlation between
#'       \code{x \%*\% s1} and \code{x \%*\% s2}, without centering;
#'       requires \code{x}.}
#'     \item{\code{"sine"}}{Sine discrepancy based on principal angles.}
#'   }
#' @param x A matrix of covariate values, required only when
#'   \code{type = "canonical"}.
#'
#' @details
#' Let \eqn{d_1}{d1} and \eqn{d_2}{d2} be the dimensions of the first and
#' second subspaces, and let \eqn{\theta_i}{theta_i} be their principal angles.
#' Then \code{type = "sine"} returns
#' \deqn{\sqrt{\sum_{i=1}^{\min(d_1,d_2)} \sin^4(\theta_i) + \max(d_1-d_2,0)}.}{
#' sqrt(sum of fourth powers of principal-angle sines + max(d1 - d2, 0)).}
#' For equal dimensions, this is the square root of the sum of fourth powers
#' of the sines. For unequal dimensions, the result depends on argument order.
#' For example, a line in a plane gives 0 with the line first and 1 with the
#' plane first. The sine discrepancy can violate the triangle inequality even
#' at equal dimensions. Use \code{type = "dist"} for a metric.
#'
#' @return A scalar giving the comparison specified by \code{type}.
#' For \code{type = "sine"}, R can return a complex scalar with zero imaginary
#' part because of roundoff in the eigenvalue calculation.
#' @examples
#' # two spaces
#' failEDR = as.matrix(cbind(c(1, 1, 0, 0, 0, 0),
#'                           c(0, 0, 1, -1, 0, 0)))
#' B = as.matrix(cbind(c(0.1, 1.1, 0, 0, 0, 0),
#'                     c(0, 0, 1.1, -0.9, 0, 0)))
#'
#' distance(failEDR, B, "dist")
#' distance(failEDR, B, "trace")
#' distance(failEDR, B, "sine")
#'
#' N=300
#' P=6
#' dataX = matrix(rnorm(N*P), N, P)
#' distance(failEDR, B, "canonical", dataX)

distance <- function(s1, s2, type = "dist", x = NULL)
{
  if (!is.matrix(s1))
    s1 = as.matrix(s1)

  if (!is.matrix(s2))
    s2 = as.matrix(s2)

  if (nrow(s1) != nrow(s2))
    stop("Dimension P of two spaces do not match.")

  if (ncol(s1) != ncol(s2))
    warning("Dimension d of two spaces do not match.")

  type = match.arg(type, c("dist", "trace", "canonical", "sine"))

  basis <- function(s) {
    peak <- apply(abs(s), 2L, max)
    if (nrow(s) < ncol(s) || any(!is.finite(peak) | peak == 0))
      stop("Each basis must have finite, linearly independent columns.")
    decomp <- qr(sweep(s, 2L, peak, "/"), tol = max(dim(s)) * .Machine$double.eps)
    if (decomp$rank != ncol(s)) stop("Each basis must have linearly independent columns.")
    qr.Q(decomp)
  }
  q1 <- basis(s1)
  q2 <- basis(s2)
  if (type == "dist")
    return(sqrt(sum((tcrossprod(q1) - tcrossprod(q2))^2)))
  if (type == "trace")
    return(sum(crossprod(q1, q2)^2) / ncol(s1))

  if (type == "canonical")
  {
    if (is.null(x))
      stop("x must be specified if use type = 'canonical'")

    if (ncol(x)!= nrow(s1))
      stop("Dimension of x is not correct.")



    return( mean(cancor(x %*% q1, x %*% q2, xcenter = FALSE, ycenter = FALSE)$cor) )
  }
  
  if (type == "sine")
  {
    Mat_1 = tcrossprod(q1)
    Mat_2 = tcrossprod(q2)
    d = eigen(Mat_1 %*% (diag(1, ncol(Mat_2)) - Mat_2))$values
    
    return( sqrt(sum(d^2)) )
  }
}
