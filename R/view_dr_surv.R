#' @title Plot Survival Data Against One Reduced Predictor
#' @name view_dr_surv
#' @export
#' @description
#' Plots observed survival times against \eqn{xB} in two dimensions, or a
#' kernel-smoothed survival surface over time and \eqn{xB} in three dimensions.
#' @param x A numeric matrix of continuous predictors, with one row per
#'   observation. Columns are not scaled internally.
#' @param y Observed follow-up times, one finite value per row of \code{x}.
#' @param censor A numeric or logical vector or one-column matrix, with one
#'   event indicator per observation (1 or `TRUE` = event, 0 or `FALSE` =
#'   censored). Missing and other values cause an error.
#' @param B A vector or one-column matrix of direction weights, with one entry
#'   per predictor. Only one reduced predictor can be displayed.
#' @param bw Kernel bandwidth for the 3D survival estimate, in units of
#'   \eqn{xB}. The default is \code{silverman(1, length(y)) * sd(x \%*\% B)}.
#' @param FUN A function, or its name, applied to \code{y} before plotting.
#'   The default logarithm requires positive times. Use \code{identity} to
#'   retain the supplied time scale. The transformation must give one finite
#'   display coordinate per observation.
#' @param type Plot type: \code{"2D"} or \code{"3D"}.
#' @param legend.add Logical; add a legend to the 2D plot.
#' @param xlab Label for the reduced-predictor axis.
#' @param ylab Label for the transformed-time axis.
#' @param zlab Label for the survival-probability axis in a 3D plot.
#' @return Invisible \code{NULL}. Called for its side-effect (producing a plot).
#' @details
#' The 3D time grid combines the observed follow-up times at empirical
#' percentiles 0 through 100 (including censored observations) with every
#' distinct event time. Percentiles use the inverse empirical distribution,
#' so each grid point is an observed time. The supplied \code{FUN} is applied
#' once to all observation times, and the corresponding transformed values
#' define the plotted grid. Risk sets use the original follow-up times.
#' When several knots share a display coordinate, the surface shows survival
#' at the latest of those follow-up times. A single distinct transformed time is displayed
#' as a line across the reduced predictor.
#'
#' @references
#' Sun, Q., Zhu, R., Wang, T. and Zeng, D. (2019). Counting process-based
#' dimension reduction methods for censored outcomes.
#' \emph{Biometrika}, 106(1), 181--196. \doi{10.1093/biomet/asy064}.
#'
#' @seealso \code{\link{orthoDr_surv}}, \code{\link{distance}}
#' @examples
#' # generate some survival data
#' N = 100; P = 4; dataX = matrix(rnorm(N*P), N, P)
#' Y = exp(-1 + dataX[,1] + rnorm(N))
#' Censor = rbinom(N, 1, 0.8)
#'
#' orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 1, method = "dm")
#' view_dr_surv(dataX, Y, Censor, orthoDr.fit$B)

view_dr_surv <- function(x, y, censor, B = NULL, bw = NULL, FUN = "log", type = "2D", legend.add = TRUE, xlab = "Reduced Direction", ylab = "Time", zlab = "Survival")
{

  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  check_survival_time(y, nrow(x))
  check_survival_censor(censor, nrow(x))

  if (is.null(B)) stop("B must be given")

  if (length(B) != ncol(x)) stop("Dimension of B does not match x")

  FUN <- match.fun(FUN)

  y2 = FUN(y)
  check_survival_time(y2, nrow(x), name = "FUN(y)")

  bx = x %*% B

  if (type == "2D")
  {
    par(mar = c(4, 4.2, 2, 2))
    plot(bx, y2, col = ifelse(censor == 1, "blue", "red"), pch = ifelse(censor == 1, 19, 3), xlab= xlab, ylab= ylab, cex.lab = 1.5)

    if (legend.add) legend("topright", c("failure", "censored"), pch = c(19, 3), col = c("blue", "red"), cex = 1.5)
  }

  if (type == "3D")
  {
    y2 = as.numeric(y2)

    if(is.null(bw))
      bw = silverman(1, length(bx))*sd(bx)

    if (!is.numeric(bw))
    {
      warning("bw must be a number")
      bw = silverman(1, length(bx))*sd(bx)
    }

    # Integer order ranks avoid rounding an empirical percentile off a jump.
    ranks = pmax(1, 1 + floor(((0:100) * as.double(length(y)) - 1) / 100))
    knots = unique(c(order(y)[ranks], which(censor == 1)))
    knots = knots[order(y2[knots], y[knots])]
    knots = knots[!duplicated(y2[knots], fromLast = TRUE)]
    timegrid = y2[knots]
    xgrid = seq(min(bx), max(bx), length.out = 100)

    log_weights = -(outer(xgrid, drop(bx), "-") / bw / sqrt(2))^2
    if (any(rowSums(is.finite(log_weights)) == 0))
      stop("Cannot compute finite Gaussian surface weights; increase bw")
    km = weighted_km(log_weights, y, censor)
    block_times = y[km$order[!duplicated(km$blocks)]]
    S = t(km$surv[match(y[knots], block_times), , drop = FALSE])

    M = mesh(xgrid, timegrid)
    colorlut <- rainbow(102, start = 0.1)

    yscale = max(y2) - min(y2)
    if (yscale == 0) yscale = 1
    xscale = max(bx) - min(bx)

    if (length(timegrid) == 1L)
      rgl::lines3d(xgrid/xscale, rep(timegrid/yscale*1.5, length(xgrid)),
                   S[, 1], col = colorlut[S[, 1]*100+1], alpha = 0.9)
    else
      surface3d(M$x/xscale, M$y/yscale*1.5, S, col = colorlut[S*100+1], alpha = 0.9, theta = 50, phi = 20, labels = c("x", "y", "z"))
    box3d(expand = 1.1, draw_front = FALSE)

    axis3d(edge = "x-+", at = seq(min(bx), max(bx), length.out = 6)/xscale,
           labels = round(seq(min(bx), max(bx), length.out = 6), 2),
           tick = TRUE, line = 0, nticks = 5, cex = 1.5, adj = c(0, 0.75))

    time_ticks = if (length(timegrid) == 1L) timegrid else
      seq(min(timegrid), max(timegrid), length.out = 6)
    axis3d(edge = "y+", at = time_ticks/yscale*1.5,
           labels = round(time_ticks, 2),
           tick = TRUE, line = 0, nticks = 6, cex = 1.5, adj = c(0, -0.25))

    axis3d(edge = "z+", at = seq(0, 1, 0.2), labels = seq(0, 1, 0.2),
           tick = TRUE, line = 1, nticks = 5, cex = 1.5, adj = 0)

    mtext3d(text=xlab, edge='x-+', line=2, cex = 1.5)
    mtext3d(text=ylab, edge='y+-', line=1.5, cex = 1.5)
    mtext3d(text=zlab, edge='z+', line=2, cex = 1.5)
  }
}
