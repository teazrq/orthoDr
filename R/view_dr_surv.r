#' @title 2D or 2D view of survival data on reduced dimension
#' @name view_dr_surv
#' @description Produce 2D or 3D plots of right censored survival data based on a given dimension reduction space
#' @param x A matrix or data.frame for features (continuous only). The algorithm will not scale the columns to unit variance
#' @param y A vector of observed time
#' @param censor A vector of censoring indicator
#' @param B The dimension reduction subspace, can only be 1 dimensional
#' @param bw A Kernel bandwidth (3D plot only) for approximating the survival function, default is the Silverman's formula
#' @param scale Using a transformed scale of the time \code{y}. Can be \code{none}, \code{log}, \code{loglog} or \code{root}. Default is \code{log}.
#' @param type \code{2D} or \code{3D} plot
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
#' view_dr_surv(dataX, Y, Censor, B = failEDR[,1], scale = "log", type = "2D")

view_dr_surv <- function(x, y, censor, B = NULL, bw = NULL, scale = c("log"), type = "2D")
{

  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y) | nrow(x) != length(censor)) stop("Number of observations do not match")

  if (is.null(B)) stop("B must be given")

  if (length(B) != ncol(x)) stop("Dimention of B does not match x")

  match.arg(scale, c("none", "log", "loglog", "root"))

  if (scale == "none") y2 = y
  if (scale == "log") y2 = log(1+y)
  if (scale == "loglog") y2 = log(1 + log(1+y))
  if (scale == "root") y2 = sqrt(y + 1)

  bx = x %*% B

  if (type == "2D")
  {
    plot.new()
    par(mar = c(4, 4.2, 2, 2))
    plot(bx, y2, col = ifelse(censor == 1, "blue", "red"), pch = ifelse(censor == 1, 19, 3), xlab = "Reduced Direction", ylab = "Time", cex.lab = 1.5)
    legend("topright", c("censored", "failure"), pch = c(19, 3), col = c("blue", "red"), cex = 1.5)
  }

  if (type == "3D")
  {

    if(is.null(bw))
      bw = silverman(1, length(bx))*sd(bx)

    if (!is.numeric(bw))
    {
      warning("bw must be a number")
      bw = silverman(1, length(bx))*sd(bx)
    }

    timegrid = seq(0, max(y2), length.out = 100)
    xgrid = seq(min(bx), max(bx), length.out = 100)

    S = matrix(NA, length(xgrid), length(timegrid))

    for (i in 1:nrow(S))
    {
      dif = xgrid[i] - bx
      k = exp(-(dif/bw)^2)
      fit = survfit(Surv(y2, censor) ~ 1, weights = k)
      S[i,] = summary(fit, times= timegrid)$surv
    }


    M = mesh(xgrid, timegrid)
    colorlut <- rainbow(102, start = 0.1)

    yscale = max(y2) - min(y2)
    xscale = max(bx) - min(bx)

    surface3d(M$x/xscale, M$y/yscale*1.5, S, col = colorlut[S*100+1], alpha = 0.9, theta = 50, phi = 20, labels = c("x", "y", "z"))
    box3d(expand = 1.1, draw_front = FALSE)

    axis3d(edge = "x-+", at = seq(min(bx), max(bx), length.out = 6)/xscale,
           labels = round(seq(min(bx), max(bx), length.out = 6), 2),
           tick = TRUE, line = 0, nticks = 5, cex = 1.5, adj = c(0, 0.75))

    axis3d(edge = "y+", at =  seq(0, max(y2) - min(y2), length.out = 6)/yscale,
           labels = round(seq(0, max(y2) - min(y2), length.out = 6), 2),
           tick = TRUE, line = 0, nticks = 6, cex = 1.5, adj = c(0, -0.25))

    axis3d(edge = "z+", at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
           tick = TRUE, line = 1, nticks = 5, cex = 1.5, adj = 0)

    mtext3d(text="Time", edge='y+-', line=1.5, cex = 1.5)
    mtext3d(text="Reduced Direction", edge='x-+', line=2, cex = 1.5)
    mtext3d(text="Survival", edge='z+', line=2, cex = 1.5)
  }
}
