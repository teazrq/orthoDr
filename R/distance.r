#' @title distance correlation
#' @description Calculate the distance correlation between two linear spaces
#' @param s1 first space
#' @param s2 second space
#' @param type type of distance measures: "dist" (default), "trace" or "canonical"
#' @param x the covairate values, for canonical correlation only
#' @return The distance between \code{s1} and \code{s2}.
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

  if (!is.vector(s2))
    s2 = as.matrix(s2)

  if (ncol(s1) != ncol(s2) || nrow(s1) != nrow(s2))
    stop("Dimention of two spaces do not match. Try using matrices.")

  match.arg(type, c("dist", "trace", "canonical"))

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
      stop("dimention of x is not correct.")

    return( mean(cancor(x %*% s1, x %*% s2, xcenter = FALSE, ycenter = FALSE)$cor) )
  }
}
