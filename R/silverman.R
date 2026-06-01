#' @title A simple Silverman bandwidth formula
#' @name silverman
#' @export
#' @description
#' Compute the Silverman rule-of-thumb bandwidth for kernel density estimation.
#'
#' @param d Integer; number of dimensions.
#' @param n Integer; number of observations.
#' @return A numeric scalar: the Silverman bandwidth.
#' @examples
#' silverman(1, 300)

silverman <- function(d, n)
{
  return((4/(d+2))^(1/(d+4))*n^(-1/(d+4)))
}
