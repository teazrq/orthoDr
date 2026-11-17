#' @title Silverman Rule-of-Thumb Bandwidth
#' @name silverman
#' @export
#' @description
#' Computes a Silverman rule-of-thumb bandwidth factor from the number of
#' dimensions and observations. It does not estimate a scale from the data.
#'
#' @param d Integer; number of dimensions.
#' @param n Integer; number of observations.
#' @return A numeric scalar, \eqn{\{4/(d+2)\}^{1/(d+4)} n^{-1/(d+4)}}.
#'   Multiply by the relevant standard deviation when a bandwidth in the
#'   original measurement units is needed.
#' @examples
#' silverman(1, 300)

silverman <- function(d, n)
{
  return((4/(d+2))^(1/(d+4))*n^(-1/(d+4)))
}
