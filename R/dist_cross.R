#' @title Gaussian Kernel Weights Between Two Data Matrices
#' @export
#' @description
#' Computes Gaussian kernel weights between every row of \code{x1} and every
#' row of \code{x2}, using their supplied coordinates without scaling.
#' @param x1 A numeric matrix with one row per observation.
#' @param x2 A numeric matrix with the same number and order of columns as \code{x1}.
#' @return A matrix with \code{nrow(x1)} rows and \code{nrow(x2)} columns.
#'   Entry \eqn{(i,j)} is \eqn{\exp(-\|x_{1i}-x_{2j}\|^2)}, where
#'   \eqn{x_{1i}} and \eqn{x_{2j}} are the corresponding rows. These are
#'   unnormalized similarity weights, with value one for identical rows.
#' @seealso \code{\link{kernel_weight}}
#' @examples
#' # two matrices
#' set.seed(1)
#' x1 = matrix(rnorm(10), 5, 2)
#' x2 = matrix(rnorm(6), 3, 2)
#' dist_cross(x1, x2)

dist_cross <- function(x1, x2)
{
  if (!is.matrix(x1) || !is.numeric(x1)) stop("x1 must be a numerical matrix")
  if (!is.matrix(x2) || !is.numeric(x2)) stop("x2 must be a numerical matrix")
  if (ncol(x1) != ncol(x2) ) stop("x1 and x2 must have the same number of columns")

  return( KernelDist_cross(x1, x2) )
}
