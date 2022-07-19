#' Counting process based sliced inverse regression model
#' 
#' The CP-SIR model for right-censored survival outcome. This model is correct
#' only under very strong assumptions, however, since it only requires an SVD,
#' the solution is used as the initial value in the orthoDr optimization.
#' 
#' @param x      A matrix for features (continuous only).
#' @param y      A vector of observed time.
#' @param censor A vector of censoring indicator.
#' @param bw     Kernel bandwidth for nonparametric estimations (one-dimensional), 
#'               the default is using Silverman's formula.
#'               
#' @return 
#' A `list` consisting of
#' 
#' \item{values}{The eigenvalues of the estimation matrix}
#' \item{vectors}{The estimated directions, ordered by eigenvalues}
#' 
#' @export
#' 
#' @references 
#' Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2019) "Counting Process Based Dimension Reduction Method for Censored Outcomes." 
#' Biometrika, 106(1), 181-196.
#' DOI: \doi{10.1093/biomet/asy064}
#' 
#' @examples
#' # This is setting 1 in Sun et. al. (2017) with reduced sample size
#' library(MASS)
#' set.seed(1)
#' N <- 200
#' P <- 6
#' V <- 0.5^abs(outer(1:P, 1:P, "-"))
#' dataX <- as.matrix(mvrnorm(N, mu = rep(0, P), Sigma = V))
#' failEDR <- as.matrix(c(1, 0.5, 0, 0, 0, rep(0, P - 5)))
#' censorEDR <- as.matrix(c(0, 0, 0, 1, 1, rep(0, P - 5)))
#' T <- rexp(N, exp(dataX %*% failEDR))
#' C <- rexp(N, exp(dataX %*% censorEDR - 1))
#' ndr <- 1
#' Y <- pmin(T, C)
#' Censor <- (T < C)
#'
#' # fit the model
#' cpsir.fit <- CP_SIR(dataX, Y, Censor)
#' distance(failEDR, cpsir.fit$vectors[, 1:ndr, drop = FALSE], "dist")
CP_SIR <- function(x, y, censor, bw = silverman(1, length(y))) {
  if (!is.matrix(x)) {
    stop("X must be a matrix")
  }
  if (nrow(x) != length(y) |
    nrow(x) != length(censor)) {
    stop("Number of observations do not match")
  }

  N <- nrow(x)
  P <- ncol(x)

  X_cov <- cov(x)

  ee <- eigen(X_cov)
  X_cov_RI <- ee$vectors %*% diag(1 / sqrt(ee$values)) %*% t(ee$vectors)

  X_sd <- scale(x, scale = FALSE) %*% X_cov_RI

  FailInd <- cbind(y, censor, "obs" = 1:N)
  FailInd <- FailInd[FailInd[, 2] == 1, ]
  OrderedFailObs <- FailInd[order(FailInd[, 1]), 3]

  timepoints <- sort(y[censor == 1])
  nFail <- length(timepoints)

  inRisk <- matrix(NA, nrow(x), nFail)

  for (i in 1:nFail) {
    inRisk[, i] <- (y >= timepoints[i])
  }

  Failure <- matrix(NA, N, nFail)
  width <- sum(censor) * bw / 2

  for (i in 1:nFail)
  {
    Failure[, i] <- (y >= timepoints[max(1, i - width)]) &
      (y <= timepoints[min(i + width, nFail)]) & (censor == 1)
  }

  # get Gu

  Gu <- matrix(NA, nFail, P)

  for (i in 1:nFail)
  {
    Gu[i, ] <- colMeans(X_sd[Failure[, i], , drop = FALSE])
  }

  Gu_risk <- matrix(NA, nFail, P)

  for (i in 1:nFail)
  {
    Gu_risk[i, ] <- colMeans(X_sd[inRisk[, i], , drop = FALSE])
  }

  Meigen <- svd(t(X_sd[OrderedFailObs, ] - Gu_risk) %*% (Gu - Gu_risk))

  return(list(
    "values" = Meigen$d,
    "vectors" = apply(X_cov_RI %*% Meigen$v, 2, function(x) {
      x / sqrt(sum(x^2))
    })
  ))
}
