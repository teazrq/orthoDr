#' Silverman's rule of thumb
#' 
#' A simple Silverman's rule of thumb bandwidth calculation.
#' 
#' @param d Number of dimension
#' @param n Number of observation
#' 
#' @return
#' A simple bandwidth choice
#' 
#' @export
#' 
#' @examples
#' silverman(1, 300)
silverman <- function(d, n) {
  return((4 / (d + 2))^(1 / (d + 4)) * n^(-1 / (d + 4)))
}
