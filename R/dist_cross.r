#' @title dist_cross
#' @description Cross distance matrix. An extension to the dist() function. Calculate the Gaussian kernel distance between rows of X1 and rows of X2
#' @param x1 first data matrix
#' @param x2 second data matrix
#' @return A distance matrix, with its (i, j)th element being the Gaussian kernel distance between ith row of \code{X1} jth row of \code{X2}.
#' @examples
#' # two matrices
#' set.seed(1)
#' x1 = matrix(rnorm(10), 5, 2)
#' x2 = matrix(rnorm(6), 3, 2)
#' dist_cross(x1, x2)

dist_cross <- function(x1, x2)
{
  if (!is.matrix(x1) | !is.numeric(x1)) stop("x1 must be a numerical matrix")
  if (!is.matrix(x2) | !is.numeric(x2)) stop("x2 must be a numerical matrix")
  if (ncol(x1) != ncol(x2) ) stop("x1 and x2 must have the same number of columns")

  return( KernelDist_cross(x1, x2) )
}
