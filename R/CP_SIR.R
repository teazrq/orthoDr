#' @title Counting-Process Based Sliced Inverse Regression
#' @name CP_SIR
#' @export
#'
#' @description
#' Estimates counting-process-based sliced inverse regression (CP-SIR)
#' directions for right-censored survival outcomes using a singular value
#' decomposition. Statistical recovery requires strong assumptions. The
#' directions are used by default to initialize \code{\link{orthoDr_surv}};
#' that function also accepts a custom starting matrix.
#'
#' @param x A finite numeric matrix of continuous predictors. Centered columns
#'   must be linearly independent. Remove constant or redundant predictors.
#' @param y Observed follow-up times, one finite value per subject, including
#'   censored subjects.
#' @param censor A numeric or logical vector or one-column matrix, with one
#'   event indicator per observation (1 or `TRUE` = event, 0 or `FALSE` =
#'   censored). Missing and other values cause an error. At least one event
#'   is required.
#' @param bw Numeric; kernel bandwidth for nonparametric estimation
#'   (one-dimensional). Default: Silverman's rule via \code{\link{silverman}}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{values}}{Singular values of the estimation matrix.}
#'   \item{\code{vectors}}{Estimated directions (columns), ordered by
#'     decreasing singular value.}
#' }
#'
#' @details
#' The estimation matrix pairs each event subject's deviation from the risk-set
#' mean with a contrast estimated from neighboring events. Smoothing the
#' neighborhood can make this matrix asymmetric. CP-SIR uses its left singular
#' vectors and maps them back to the original feature coordinates.
#'
#' Events at the same observed time share the midpoint of their event ranks
#' when selecting a neighborhood. Neighborhood endpoints are bounded by the
#' first and last event ranks and truncated to integer indices. Risk sets
#' include all observations at or after the event time.
#'
#' @references
#' Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2019). Counting process-based
#' dimension reduction methods for censored outcomes.
#' \emph{Biometrika}, 106(1), 181--196. \doi{10.1093/biomet/asy064}.
#'
#' @seealso \code{\link{orthoDr_surv}}, \code{\link{silverman}}
#' @examples
#' # Setting 1 from Sun et al. (2019), reduced sample size
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
  check_survival_time(y, nrow(x))
  check_survival_censor(censor, nrow(x))

  N = nrow(x)
  P = ncol(x)
  nFail = sum(censor == 1)
  if (nFail == 0) stop("CP_SIR requires at least one event")

  coords = survival_coordinates(x)
  design = survival_design_svd(coords$x, nv = P)
  whitener = tcrossprod(sweep(design$v, 2L, sqrt(N - 1L) / design$d, "*"), design$v)
  X_sd = coords$x %*% whitener

  events = which(censor == 1)
  OrderedFailObs = events[order(y[events])]
  timepoints = y[OrderedFailObs]
  nFail = length(timepoints)

  inRisk = matrix(NA, nrow(x), nFail)

  for (i in 1:nFail)
    inRisk[, i] = (y >= timepoints[i])

  Failure = matrix(NA, N, nFail)
  width = sum(censor)*bw/2
  # Tied event times share their midpoint rank. Clamp before positive-index
  # truncation to retain the original neighborhood rule for untied times.
  event.rank = rank(timepoints, ties.method = "average")

  for (i in 1:nFail)
  {
    lower = as.integer(max(1, event.rank[i] - width))
    upper = as.integer(min(nFail, event.rank[i] + width))
    Failure[, i] = (y >= timepoints[lower]) & (y <= timepoints[upper]) & (censor == 1)
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

  Meigen = svd(t(X_sd[OrderedFailObs, , drop = FALSE] - Gu_risk) %*% (Gu - Gu_risk))
  # The event subject's residual is on the left of the moment.
  vectors = raw_predictor_directions(whitener %*% Meigen$u, coords$scale)

  return(list("values" = Meigen$d,
              "vectors" = vectors))
}
