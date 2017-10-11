#' @title Counting process based sliced inverse regression model
#' @name CP_SIR
#' @description The CP-SIR model for right-censored survival outcome. This model is correct only under very strong assumptions, however, since it only requires an SVD, the solution is used as the initial value in the orthoDr optimization.
#' @param x A matrix for features (continuous only).
#' @param y A vector of observed time.
#' @param censor A vector of censoring indicator.
#' @param bw Kernel bandwidth for nonparametric estimations (one-dimensional), the default is using Silverman's formula.
#' @return A list consisting of
#' \item{values}{The eigenvalues of the estimation matrix}
#' \item{vectors}{The estimated directions, ordered by eigenvalues}
#' @references Sun, Q., Zhu, R., Wang T. and Zeng D. "Counting Process Based Dimension Reduction Method for Censored Outcomes." (2017)
#' \url{https://arxiv.org/abs/1704.05046} .
#' @examples
#' # generate some survival data
#' set.seed(1)
#' N = 350
#' P = 6
#' dataX = matrix(rnorm(N*P), N, P)
#' failEDR = as.matrix(cbind(c(1, 1, 0, 0, 0, 0, rep(0, P-6)),
#'                           c(0, 0, 1, -1, 0, 0, rep(0, P-6))))
#' censorEDR = as.matrix(c(0, 1, 0, 1, 1, 1, rep(0, P-6)))
#' T = exp(-2.5 + dataX %*% failEDR[,1] +
#'     0.5*(dataX %*% failEDR[,1])*(dataX %*% failEDR[,2]) + 0.25*log(-log(1-runif(N))))
#' C = exp( -0.5 + dataX %*% censorEDR  + log(-log(1-runif(N))))
#' Y = pmin(T, C)
#' Censor = (T < C)
#'
#' # fit the model
#' CP_SIR.fit = CP_SIR(dataX, Y, Censor)
#'
#' # compare with the true direction
#' distance(failEDR, CP_SIR.fit$vectors[, 1:2], "dist")

CP_SIR <- function(x, y, censor, bw = silverman(1, length(y)))
{
  if (!is.matrix(x)) stop("X must be a matrix")
  if (nrow(x) != length(y) | nrow(x) != length(censor)) stop("Number of observations do not match")

  N = nrow(x)
  P = ncol(x)

  X_cov = cov(x)

  ee = eigen(X_cov)
  X_cov_RI = ee$vectors %*% diag(1/sqrt(ee$values)) %*% t(ee$vectors)
  # X_cov_R = ee$vectors %*% diag(sqrt(ee$values)) %*% t(ee$vectors)

  X_sd = scale(x, scale = FALSE) %*% X_cov_RI

  FailInd = cbind(y, censor, "obs" = 1:N)
  FailInd = FailInd[FailInd[, 2] == 1, ]
  OrderedFailObs = FailInd[order(FailInd[,1]), 3]

  timepoints = sort(y[censor == 1])
  nFail = length(timepoints)

  inRisk = matrix(NA, nrow(x), nFail)

  for (i in 1:nFail)
    inRisk[, i] = (y >= timepoints[i])


  Failure = matrix(NA, N, nFail)
  width = sum(censor)*bw/2

  for (i in 1:nFail)
  {
    Failure[, i] = (y >= timepoints[max(1, i-width)]) & (y <= timepoints[min(i+width, nFail)]) & (censor == 1)
  }

  # get Gu

  Gu = matrix(NA, nFail, P)

  for (i in 1:nFail)
  {
    Gu[i, ] = colMeans(X_sd[Failure[, i], , drop = FALSE ])
  }

  Gu_risk = matrix(NA, nFail, P)

  for (i in 1:nFail)
  {
    Gu_risk[i, ] = colMeans(X_sd[inRisk[, i], , drop = FALSE ])
  }

  Meigen = svd(t(X_sd[OrderedFailObs, ] - Gu_risk) %*% (Gu - Gu_risk))

  return(list("values" = Meigen$d,
              "vectors" = apply(X_cov_RI %*% Meigen$v, 2, function(x) {x/sqrt(sum(x^2))})))
}
