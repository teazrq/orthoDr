#' 2D or 2D view of survival data on reduced dimension
#' 
#' Produce 2D or 3D plots of right censored survival data based on a given
#' dimension reduction space
#' 
#' @param x          A `matrix` or `data.frame` for features (continuous only).
#'                   The algorithm will not scale the columns to unit variance
#' @param y          A `vector` of observed time
#' @param censor     A `vector` of censoring indicator
#' @param B          The dimension reduction subspace, can only be 1 dimensional
#' @param bw         A Kernel bandwidth (3D plot only) for approximating the 
#'                   survival function, default is the Silverman's formula
#' @param FUN        A scaling function applied to the time points `y`.
#'                   Default is `"log"`.
#' @param type       `2D` or `3D` plot
#' @param legend.add Should legend be added (2D plot only)
#' @param xlab       x axis label
#' @param ylab       y axis label
#' @param zlab       z axis label
#' 
#' @return 
#' 
#' An `rgl` object that is rendered.
#' 
#' @export
#' 
#' @importFrom grDevices rainbow
#' @importFrom graphics legend par plot plot.new
#' @importFrom rgl axis3d mtext3d box3d surface3d
#' @importFrom plot3D mesh
#' 
#' @references 
#' Sun, Q., Zhu, R., Wang, T., & Zeng, D. (2019). Counting process-based dimension reduction methods for censored outcomes. Biometrika, 106(1), 181-196.
#' DOI: \doi{10.1093/biomet/asy064}
#' 
#' @examples
#' # generate some survival data
#' N <- 100
#' P <- 4
#' dataX <- matrix(rnorm(N * P), N, P)
#' Y <- exp(-1 + dataX[, 1] + rnorm(N))
#' Censor <- rbinom(N, 1, 0.8)
#'
#' orthoDr.fit <- orthoDr_surv(dataX, Y, Censor, ndr = 1, method = "dm")
#' view_dr_surv(dataX, Y, Censor, orthoDr.fit$B)
view_dr_surv <- function(x, y, censor, B = NULL, bw = NULL, FUN = "log", type = "2D", legend.add = TRUE, xlab = "Reduced Direction", ylab = "Time", zlab = "Survival") {
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y) | nrow(x) != length(censor)) stop("Number of observations do not match")

  if (is.null(B)) stop("B must be given")

  if (length(B) != ncol(x)) stop("Dimension of B does not match x")

  FUN <- match.fun(FUN)

  y2 <- FUN(y)

  bx <- x %*% B

  if (type == "2D") {
    plot.new()
    par(mar = c(4, 4.2, 2, 2))
    plot(bx, y2, col = ifelse(censor == 1, "blue", "red"), pch = ifelse(censor == 1, 19, 3), xlab = xlab, ylab = ylab, cex.lab = 1.5)

    if (legend.add) legend("topright", c("failure", "censored"), pch = c(19, 3), col = c("blue", "red"), cex = 1.5)
  }

  if (type == "3D") {
    if (is.null(bw)) {
      bw <- silverman(1, length(bx)) * sd(bx)
    }

    if (!is.numeric(bw)) {
      warning("bw must be a number")
      bw <- silverman(1, length(bx)) * sd(bx)
    }

    timegrid <- sort(unique(c(y2, seq(0, max(y2), length.out = 100))))
    xgrid <- seq(min(bx), max(bx), length.out = 100)

    S <- matrix(NA, length(xgrid), length(timegrid))

    for (i in 1:nrow(S))
    {
      dif <- xgrid[i] - bx
      k <- exp(-0.5 * (dif / bw)^2)
      fit <- survfit(Surv(y2, censor) ~ 1, weights = k)
      S[i, ] <- summary(fit, times = timegrid)$surv
    }


    M <- mesh(xgrid, timegrid)
    colorlut <- rainbow(102, start = 0.1)

    yscale <- max(y2) - min(y2)
    xscale <- max(bx) - min(bx)

    surface3d(M$x / xscale, M$y / yscale * 1.5, S, col = colorlut[S * 100 + 1], alpha = 0.9, theta = 50, phi = 20, labels = c("x", "y", "z"))
    box3d(expand = 1.1, draw_front = FALSE)

    axis3d(
      edge = "x-+", at = seq(min(bx), max(bx), length.out = 6) / xscale,
      labels = round(seq(min(bx), max(bx), length.out = 6), 2),
      tick = TRUE, line = 0, nticks = 5, cex = 1.5, adj = c(0, 0.75)
    )

    axis3d(
      edge = "y+", at = seq(0, max(y2) - min(y2), length.out = 6) / yscale,
      labels = round(seq(0, max(y2) - min(y2), length.out = 6), 2),
      tick = TRUE, line = 0, nticks = 6, cex = 1.5, adj = c(0, -0.25)
    )

    axis3d(
      edge = "z+", at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
      tick = TRUE, line = 1, nticks = 5, cex = 1.5, adj = 0
    )

    mtext3d(text = xlab, edge = "x-+", line = 2, cex = 1.5)
    mtext3d(text = ylab, edge = "y+-", line = 1.5, cex = 1.5)
    mtext3d(text = zlab, edge = "z+", line = 2, cex = 1.5)
  }
}
