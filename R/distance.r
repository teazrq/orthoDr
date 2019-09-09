#' Compute Distance Correlation
#'
#' Calculate the distance correlation between two linear spaces.
#' 
#' @param s1   First space
#' @param s2   Second space
#' @param type Type of distance measures: `"dist"` (default), `"trace"`,
#'             `"canonical"` or `"sine"`
#' @param x    The covariate values, for canonical correlation only.
#' 
#' @return 
#' The distance between `s1` and `s2`.
#' 
#' @export 
#' 
#' @examples
#' # two spaces
#' failEDR <- as.matrix(cbind(
#'   c(1, 1, 0, 0, 0, 0),
#'   c(0, 0, 1, -1, 0, 0)
#' ))
#' B <- as.matrix(cbind(
#'   c(0.1, 1.1, 0, 0, 0, 0),
#'   c(0, 0, 1.1, -0.9, 0, 0)
#' ))
#'
#' distance(failEDR, B, "dist")
#' distance(failEDR, B, "trace")
#'
#' N <- 300
#' P <- 6
#' dataX <- matrix(rnorm(N * P), N, P)
#' distance(failEDR, B, "canonical", dataX)
distance <- function(s1, s2, type = "dist", x = NULL) {
  if (!is.matrix(s1)) {
    s1 <- as.matrix(s1)
  }

  if (!is.vector(s2)) {
    s2 <- as.matrix(s2)
  }

  if (nrow(s1) != nrow(s2)) {
    stop("Dimension P of two spaces do not match.")
  }

  if (ncol(s1) != ncol(s2)) {
    warning("Dimention d of two spaces do not match.")
  }

  match.arg(type, c("dist", "trace", "canonical", "sine"))

  if (type == "dist") {
    Mat_1 <- s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 <- s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    return(sqrt(sum((Mat_1 - Mat_2)^2)))
  }

  if (type == "trace") {
    Mat_1 <- s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 <- s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    return(sum(diag(x = Mat_1 %*% Mat_2)) / ncol(s1))
  }

  if (type == "canonical") {
    if (is.null(x)) {
      stop("x must be specified if use type = 'canonical'")
    }

    if (ncol(x) != nrow(s1)) {
      stop("Dimension of x is not correct.")
    }

    return(mean(cancor(x %*% s1, x %*% s2, xcenter = FALSE, ycenter = FALSE)$cor))
  }

  if (type == "sine") {
    Mat_1 <- s1 %*% solve(t(s1) %*% s1) %*% t(s1)
    Mat_2 <- s2 %*% solve(t(s2) %*% s2) %*% t(s2)
    d <- eigen(Mat_1 %*% (diag(1, ncol(Mat_2)) - Mat_2))$values

    return(sqrt(sum(d^2)))
  }
}
