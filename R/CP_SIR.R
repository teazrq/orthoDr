#' @title Counting-Process Based Sliced Inverse Regression
#' @name CP_SIR
#' @export
#'
#' @description
#' Fits the CP-SIR model for right-censored survival outcomes. This model is
#' correct only under strong assumptions, but since it requires only an SVD,
#' its solution is commonly used as the initial value in
#' \code{\link{orthoDr_surv}} optimization.
#'
#' @param x A matrix (n x p) of features (continuous only).
#' @param y A numeric vector of observed survival times.
#' @param censor A numeric vector of censoring indicators (1 = event,
#'   0 = censored).
#' @param bw Numeric; kernel bandwidth for nonparametric estimation
#'   (one-dimensional). Default: Silverman's rule via \code{\link{silverman}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{values}}{Eigenvalues of the estimation matrix.}
#'   \item{\code{vectors}}{Estimated directions (columns), ordered by
#'     eigenvalues.}
#' }
#'
#' @references
#' Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2017). Counting process based
#' dimension reduction method for censored outcomes.
#' \emph{arXiv preprint} \doi{10.48550/arXiv.1704.05046}.
#'
#' @seealso \code{\link{orthoDr_surv}}, \code{\link{silverman}}
#' @examples
#' # This is setting 1 in Sun et. al. (2017) with reduced sample size
#' library(MASS)
#' set.seed(1)
#' N = 200; P = 6
#' V=0.5^abs(outer(1:P, 1:P, "-"))
#' dataX = as.matrix(mvrnorm(N, mu=rep(0,P), Sigma=V))
#' failEDR = as.matrix(c(1, 0.5, 0, 0, 0, rep(0, P-5)))
#' censorEDR = as.matrix(c(0, 0, 0, 1, 1, rep(0, P-5)))
#' T = rexp(N, exp(dataX %*% failEDR))
#' C = rexp(N, exp(dataX %*% censorEDR - 1))
#' ndr = 1
#' Y = pmin(T, C)
#' Censor = (T < C)
#'
#' # fit the model
#' cpsir.fit = CP_SIR(dataX, Y, Censor)
#' distance(failEDR, cpsir.fit$vectors[, 1:ndr, drop = FALSE], "dist")

CP_SIR <- function(x, y, censor, bw = silverman(1, length(y)))
{
  if (!is.matrix(x)) stop("X must be a matrix")
  if (nrow(x) != length(y) || nrow(x) != length(censor)) stop("Number of observations do not match")

  N = nrow(x)
  P = ncol(x)

  X_cov = cov(x)

  ee = eigen(X_cov)
  X_cov_RI = ee$vectors %*% diag(1/sqrt(ee$values)) %*% t(ee$vectors)

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
