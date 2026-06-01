#' @title Partial Sliced Averaged Variance Estimation
#' @name pSAVE
#' @export
#'
#' @description
#' Fits the partial SAVE model for dose-response settings. This model relies
#' on strong assumptions; its solution is primarily used as an initial value
#' in \code{\link{orthoDr_pdose}} optimization.
#'
#' @param x A matrix (n x p) of features (continuous only).
#' @param a A numeric vector of observed dose levels.
#' @param r A numeric vector of rewards (outcomes).
#' @param ndr Integer; number of structural dimensions to estimate.
#'   Default: \code{2}.
#' @param nslices0 Integer; number of slices used for SAVE.
#'   Default: \code{2}.
#'
#' @return A matrix whose columns are the estimated basis vectors of the
#'   central subspace, ordered by eigenvalues.
#'
#' @references
#' Feng, Z., Wen, M.X., Yu, Z. and Zhu, L. (2013). On partial sufficient
#' dimension reduction with applications to partially linear multi-index
#' models. \emph{Journal of the American Statistical Association},
#' 108(501), 237--256. \doi{10.1080/01621459.2013.849167}.
#'
#' @seealso \code{\link{orthoDr_pdose}}
#'
#' @examples
#' set.seed(1)
#' N <- 200; P <- 4
#' X <- matrix(rnorm(N * P), N, P)
#' dose <- runif(N, 0, 2)
#' reward <- X[, 1] + rnorm(N)
#' pSAVE(X, dose, reward, ndr = 1)

pSAVE <- function(x, a, r, ndr = 2, nslices0 =2){

  if (!is.matrix(x)) stop("X must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) || nrow(x) != length(a)) stop("Number of observations do not match")

  train = list(x=x,a=a,r=r)
  n = nrow(x)
  p = ncol(x)
  newtrain = train
  a = train$a
  a = sort(a)

  Z = a[(n/2-50):(n/2+50)]
  M_i = list()
  jk = 0

  for (i in Z){

    jk = jk +1
    newZ = c()
    newZ[which(train$a <= i)] = 1
    newZ[which(train$a > i)] = 0

    newtrain$a = newZ
    dimdr = dr(formula = r ~ x , data = newtrain ,group = ~a,
               nslices = nslices0, chi2approx = "bx",
               numdir = p, method = "save")
    M_i[[jk]] = dimdr$M

  }

  M_total = matrix(0,p,p)
  for (j in 1:length(Z)){
    M_total =  M_total + M_i[[j]]
  }

  Beta = eigen(M_total, symmetric =F, only.values = FALSE, EISPACK = FALSE)
  B = as.matrix(Beta$vectors[,1:ndr])
  return(B)
}


