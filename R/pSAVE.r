#' Partial Sliced Averaged Variance Estimation
#' 
#' The partial-SAVE model. This model is correct only under very strong 
#' assumptions, the solution is used as the initial value in the `orthoDr`
#' optimization.
#' 
#' @param x        A `matrix` for features (continuous only).
#' @param a        A `vector` of observed dose levels (continuous only).
#' @param r        A `vector` of reward (outcome).
#' @param ndr      The dimension structure
#' @param nslices0 Number of slides used for save
#' 
#' @return 
#' 
#' A `list` consisting of:
#' 
#' \item{vectors}{The basis of central subspace, ordered by eigenvalues}
#' 
#' @export
#' 
#' @references 
#' Feng, Z., Wen, X. M., Yu, Z., & Zhu, L. (2013). On partial sufficient dimension reduction with applications to partially linear multi-index models. 
#' Journal of the American Statistical Association, 108(501), 237-246.
#' DOI: \doi{10.1080/01621459.2012.746065}
#' 
pSAVE <- function(x, a, r, ndr = 2, nslices0 = 2) {
  if (!is.matrix(x)) stop("X must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) | nrow(x) != length(a)) stop("Number of observations do not match")

  train <- list(x = x, a = a, r = r)
  n <- nrow(x)
  p <- ncol(x)
  newtrain <- train
  a <- train$a
  a <- sort(a)

  Z <- a[(n / 2 - 50):(n / 2 + 50)]
  M_i <- list()
  jk <- 0

  for (i in Z) {
    jk <- jk + 1
    newZ <- c()
    newZ[which(train$a <= i)] <- 1
    newZ[which(train$a > i)] <- 0

    newtrain$a <- newZ
    dimdr <- dr(
      formula = r ~ x, data = newtrain, group = ~a,
      nslices = nslices0, chi2approx = "bx",
      numdir = p, method = "save"
    )
    M_i[[jk]] <- dimdr$M
  }

  M_total <- matrix(0, p, p)
  for (j in 1:length(Z)) {
    M_total <- M_total + M_i[[j]]
  }

  Beta <- eigen(M_total, symmetric = F, only.values = FALSE, EISPACK = FALSE)
  B <- as.matrix(Beta$vectors[, 1:ndr])
  return(B)
}
