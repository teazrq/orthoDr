#' Hazard Mave for Censored Survival Data
#' 
#' This is an almost direct R translation of Xia, Zhang & Xu's (2010) hMave 
#' `MATLAB` code. We implemented further options for setting a different initial
#' value. The computational algorithm does not utilize the orthogonality
#' constrained optimization.
#' 
#' @param x A matrix for features.
#' @param y A vector of observed time.
#' @param censor A vector of censoring indicator.
#' @param m0 number of dimensions to use
#' @param B0 initial value of B. This is a feature we implemented.
#' 
#' @return 
#' 
#' A `list` consisting of
#' 
#' \item{B}{The estimated B matrix}
#' \item{cv}{Leave one out cross-validation error}
#' 
#' @export
#' 
#' @references
#' Xia, Y., Zhang, D., & Xu, J. (2010). Dimension reduction and semiparametric estimation of survival models. 
#' Journal of the American Statistical Association, 105(489), 278-290.
#' DOI: \doi{10.1198/jasa.2009.tm09372}
#' 
#' @examples
#' # generate some survival data
#' set.seed(1)
#' P <- 7
#' N <- 150
#' dataX <- matrix(runif(N * P), N, P)
#' failEDR <- as.matrix(cbind(c(1, 1.3, -1.3, 1, -0.5, 0.5, -0.5, rep(0, P - 7))))
#' T <- exp(dataX %*% failEDR + rnorm(N))
#' C <- runif(N, 0, 15)
#' Y <- pmin(T, C)
#' Censor <- (T < C)
#'
#' # fit the model
#' hMave.fit <- hMave(dataX, Y, Censor, 1)
hMave <- function(x, y, censor, m0, B0 = NULL) {
  # this function returns B and cv

  # Input : x (n x p); y (n x 1); m0 (integer)
  # Output: B (p x m0)
  # m0: dimension e.g. it is 1 for single-index model

  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.matrix(y)) stop("y must be a matrix")
  if (!is.matrix(censor)) stop("censor must be a matrix")

  censor <- as.matrix(as.numeric(censor))
  if (!all(censor %in% c(0, 1))) {
    stop("censor can only have two levels 0/1")
  }

  n <- nrow(x)
  p <- ncol(x)

  onen <- matrix(1, n, 1)
  x <- scale(x, scale = FALSE)

  xeigen <- eigen(t(x) %*% x / n)

  ss <- xeigen$vectors[, rev(1:p), drop = FALSE]
  D <- diag(rev(xeigen$values))

  ss <- (ss / repmat(matrix(sqrt(diag(D)), 1, p), p, 1))

  x <- x %*% ss

  # if initial value is given
  if (!is.null(B0)) {
    if (!is.matrix(B0)) stop("B0 must be a matrix")
    if (nrow(B0) != ncol(x)) stop("Dimension of B0 is not correct")
    B0 <- solve(ss) %*% B0
    B <- B0
    m <- ncol(B)
  } else {
    m <- p
    B <- diag(p)
  }

  Ba <- B
  BI <- Ba
  B <- B[, 1:m, drop = FALSE]

  noip0 <- 0
  noip1 <- 1
  iterstop <- 0
  Btmp <- B

  rige <- sd(y) * mean(apply(x, 2, sd))
  y <- scale(y)

  # x = x./repmat(std(x,[],1),n,1);

  I <- order(y)
  y <- y[I, , drop = FALSE]
  x <- x[I, , drop = FALSE]
  censor <- censor[I, , drop = FALSE]
  censor0 <- censor

  niter <- floor(p * 3 / 2)
  # ch = (sqrt(p)/n^(1/(p+4))*n^(2/(m0+4)))^(1/niter);
  yc <- y

  ky2 <- repmat(yc, 1, n)
  ky2 <- (ky2 - t(ky2))^2
  censor <- repmat(censor, 1, n)


  for (iter in 1:niter)
  {
    adj <- p^(1 / iter) * n^(1 / (m0 + 4) - 1 / (m + 4))
    hy <- sd(yc) / n^(1 / 5) * p^(1 / (iter + 1))
    xB <- x %*% Ba
    h <- p^(1 / (iter + 1)) * mean(apply(xB, 2, sd)) / n^(1 / (m + 4))
    h2 <- 2 * h * h * adj

    H <- matrix(1, 1, n)

    for (i in 1:n)
    {
      xBi <- (xB - repmat(xB[i, , drop = FALSE], n, 1))
      k <- exp(-rowSums(xBi^2) / h2)
      kH <- pnorm((y[i] - y) / hy)
      H[, i] <- t(k) %*% kH / sum(k)
    }

    ky <- censor * exp(-ky2 / (2 * hy * hy)) / hy / repmat(H, n, 1)
    n1 <- ncol(ky)

    ABI <- matrix(0, m, n * n)

    for (iiter in 1:max(1, (m < p) * floor(m / 2)))
    {
      dd <- matrix(0, m * p, m * p)
      dc <- matrix(0, m * p, 1)

      for (j in 1:n)
      {
        xij <- x - repmat(x[j, , drop = FALSE], n, 1)
        sxij <- rowSums((xij %*% Ba)^2) + rowSums(xij^2) / 1.5^iter
        ker <- exp(-sxij / h2)
        rker <- repmat(as.matrix(ker), 1, p + 1)
        onexi <- cbind(xij %*% B, onen)
        xk <- t(onexi * rker[, 1:(m + 1)])
        abi <- solve(xk %*% onexi + diag(1, m + 1) / n) %*% (xk %*% ky)

        kxij <- t(xij * rker[, 1:p, drop = FALSE])
        kxijy <- kxij %*% (ky - repmat(abi[m + 1, , drop = FALSE], n, 1))
        ddx <- kxij %*% xij

        for (k1 in 1:m)
        {
          ka <- ((k1 - 1) * p + 1):(k1 * p)
          dc[ka, ] <- dc[ka, ] + kxijy %*% t(abi[k1, , drop = FALSE])
        }

        tmp <- abi[1:m, , drop = FALSE] %*% t(abi[1:m, , drop = FALSE])

        dd <- dd + kronecker(tmp, ddx)

        ABI[, ((j - 1) * n1 + 1):(j * n1)] <- abi[1:m, , drop = FALSE]
      }

      B <- solve(dd + rige * diag(1, length(dc)) / n) %*% dc

      B0 <- matrix(B, p, m) %*% ABI

      B <- eigen(B0 %*% t(B0))$vectors[, rev(1:p), drop = FALSE]

      B <- B[, (p - m + 1):p, drop = FALSE]
      Ba <- B

      if (max(svd(B %*% t(B) - BI %*% t(BI))$d) < 0.001) {
        break
      }

      BI <- B
    }

    mb <- m
    ma <- max(m - 1, m0)
    m <- ma
    B <- B[, (mb - m + 1):mb, drop = FALSE]
    Ba <- B

    if (max(svd(B %*% t(B) - BI %*% t(BI))$d) < 0.001 & iter > p + 3) {
      break
    }

    BI <- Ba

    if (iter == 1) {
      B00 <- B[, (ncol(B) - m0 + 1):ncol(B), drop = FALSE]
    }
  }

  h2 <- 2 * n^(-2 / (m0 + 3))

  x <- x[censor0 == 1, , drop = FALSE]
  y <- y[censor0 == 1, , drop = FALSE]
  n <- nrow(y)
  ky <- scale(y)
  ky <- exp(-(repmat(ky, 1, n) - repmat(t(ky), n, 1))^2 / h2)

  kye <- ky

  cv <- 0

  for (j in 1:n)
  {
    xij <- x - repmat(x[j, , drop = FALSE], n, 1)
    sxij <- rowSums((xij %*% Ba)^2)
    ker <- as.matrix(exp(-sxij / h2))
    ker[j] <- 0
    kye[j, ] <- t(ker) %*% ky / sum(ker)
  }

  cv <- sum((ky - kye)^2)
  B <- ss %*% B

  apply(B, 2, function(x) x / sqrt(sum(x^2)))

  return(list("B" = B, "cv" = cv))
}
