#' @title Distance Between Two Linear Subspaces
#' @name distance
#' @export
#'
#' @description
#' Compares the linear subspaces spanned by the columns of \code{s1} and
#' \code{s2} using projection-matrix distances or correlations.
#'
#' @param s1 A matrix with linearly independent columns spanning the first subspace.
#' @param s2 A matrix with linearly independent columns spanning the second subspace,
#'   with the same number of rows as \code{s1}.
#' @param type Character; the distance measure to use:
#'   \describe{
#'     \item{\code{"dist"}}{Frobenius norm of the difference between
#'       projection matrices (default).}
#'     \item{\code{"trace"}}{Trace correlation. For subspaces of equal
#'       dimension, 1 indicates identical subspaces.}
#'     \item{\code{"canonical"}}{Mean canonical correlation between
#'       \code{x \%*\% s1} and \code{x \%*\% s2}, without centering;
#'       requires \code{x}.}
#'     \item{\code{"sine"}}{Sine metric between subspaces.}
#'   }
#' @param x A matrix of covariate values, required only when
#'   \code{type = "canonical"}.
#'
#' @return A numeric scalar: the distance (or correlation) between the two
#'   subspaces. Interpretation depends on \code{type}.
#' @examples
#' # two spaces
#' failEDR = as.matrix(cbind(c(1, 1, 0, 0, 0, 0),
#'                           c(0, 0, 1, -1, 0, 0)))
#' B = as.matrix(cbind(c(0.1, 1.1, 0, 0, 0, 0),
#'                     c(0, 0, 1.1, -0.9, 0, 0)))
#'
#' distance(failEDR, B, "dist")
#' distance(failEDR, B, "trace")
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

  if (type == "dist")
  {
    Mat_1 = s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 = s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    return( sqrt(sum((Mat_1 - Mat_2)^2)) )
  }

  if (type == "trace")
  {
    Mat_1 = s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 = s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    return( sum(diag(x = Mat_1 %*% Mat_2 ))/ncol(s1) )
  }

  if (type == "canonical")
  {
    if (is.null(x))
      stop("x must be specified if use type = 'canonical'")

    if (ncol(x)!= nrow(s1))
      stop("Dimension of x is not correct.")



    return( mean(cancor(x %*% s1, x %*% s2, xcenter = FALSE, ycenter = FALSE)$cor) )
  }
  
  if (type == "sine")
  {
    Mat_1 = s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 = s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    d = eigen(Mat_1 %*% (diag(1, ncol(Mat_2)) - Mat_2))$values
    
    return( sqrt(sum(d^2)) )
  }
}
