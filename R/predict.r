#' Predictions under `orthoDr` models
#' 
#' The prediction function for `orthoDr` fitted models
#' 
#' @param object A fitted `orthoDr` object
#' @param testx  Testing data
#' @param ...    Additional parameters, not used.
#' 
#' @return 
#' The predicted object
#' 
#' @export
#' 
#' @examples
#' # generate some survival data
#' N <- 100
#' P <- 4
#' dataX <- matrix(rnorm(N * P), N, P)
#' Y <- exp(-1 + dataX[, 1] + rnorm(N))
#' Censor <- rbinom(N, 1, 0.8)
#'
#' # fit the model with keep.data = TRUE
#' orthoDr.fit <- orthoDr_surv(dataX, Y, Censor,
#'   ndr = 1,
#'   method = "dm", keep.data = TRUE
#' )
#'
#' # predict 10 new observations
#' predict(orthoDr.fit, matrix(rnorm(10 * P), 10, P))
#'
#' # generate some personalized dose scenario
#'
#' exampleset <- function(size, ncov) {
#'   X <- matrix(runif(size * ncov, -1, 1), ncol = ncov)
#'   A <- runif(size, 0, 2)
#'
#'   Edr <- as.matrix(c(0.5, -0.5))
#'
#'   D_opt <- X %*% Edr + 1
#'
#'   mu <- 2 + 0.5 * (X %*% Edr) - 7 * abs(D_opt - A)
#'
#'   R <- rnorm(length(mu), mu, 1)
#'
#'   R <- R - min(R)
#'
#'   datainfo <- list(X = X, A = A, R = R, D_opt = D_opt, mu = mu)
#'   return(datainfo)
#' }
#'
#' # generate data
#'
#' set.seed(123)
#' n <- 150
#' p <- 2
#' ndr <- 1
#' train <- exampleset(n, p)
#' test <- exampleset(500, p)
#'
#' # the direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = 0.1,
#'   method = "direct", K = as.integer(sqrt(n)), keep.data = TRUE,
#'   maxitr = 150, verbose = FALSE, ncore = 2
#' )
#'
#' predict(orthofit, test$X)
#'
#' # the pseudo direct learning method
#' orthofit <- orthoDr_pdose(train$X, train$A, train$R,
#'   ndr = ndr, lambda = seq(0.1, 0.2, 0.01),
#'   method = "pseudo_direct", K = as.integer(sqrt(n)), keep.data = TRUE,
#'   maxitr = 150, verbose = FALSE, ncore = 2
#' )
#'
#' predict(orthofit, test$X)
predict.orthoDr <- function(object, testx, ...) {
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.matrix(testx) || !is.numeric(testx)) stop("testx must be a numerical matrix")

  if (class(object)[2] != "fit") {
    stop("This is not an orthoDr fitted object")
  }

  if (!object$keep.data) {
    stop("Need the original data for prediction. Please specify keep.data = TRUE in model fitting.")
  }

  # predict survival functions on the testing data

  if (class(object)[3] == "surv") {
    pred <- predict_orthoDr_surv(object, testx, ...)
  }

  # predict regression outcome on the testing data

  if (class(object)[3] == "reg") {
    pred <- predict_orthoDr_reg(object, testx, ...)
  }

  # predict rewards and optimal dose on the testing data

  if (class(object)[3] == "pdose") {
    pred <- predict_orthoDr_pdose(object, testx, ...)
  }

  class(pred) <- c("orthoDr", "predict", class(object)[3])

  return(pred)
}

#' Internal prediction function for survival models
#' 
#' Provides support for predicting values underneath a survival model.
#' 
#' @param object Fitted object
#' @param testx  Testing data
#' @param ...    Additional parameters, not used.
#' 
#' @return 
#' The predicted object
#' 
#' @noRd
predict_orthoDr_surv <- function(object, testx, ...) {
  # transform the covariates into the same scale
  x <- object$x
  xscale <- apply(x, 2, sd)
  xmean <- apply(x, 2, mean)
  X <- scale(x)

  testX <- sweep(testx, 2, xmean, FUN = "-")
  testX <- sweep(testX, 2, xscale, FUN = "/")

  XB <- X %*% object$B
  testXB <- testX %*% object$B

  XBscale <- apply(XB, 2, sd)
  testXB <- sweep(testXB, 2, XBscale, FUN = "/")

  XB <- scale(XB) / object$bw / sqrt(2)
  testXB <- testXB / object$bw / sqrt(2)

  # kernel matrix between x and testx
  testKernel <- KernelDist_cross(testXB, XB)

  # this method does not deal with ties. I need to fix this later on

  testKernel <- testKernel[, order(object$y), drop = FALSE]
  Censor <- (object$censor[order(object$y)] == 1)

  inrisk <- apply(testKernel, 1, cumsum)
  totalweights <- inrisk[nrow(X), , drop = FALSE]
  inrisk <- sweep(-inrisk, 2, totalweights, FUN = "+")

  testKernel <- sweep(testKernel, 2, Censor, FUN = "*")
  testKernel <- t(testKernel)

  lambda <- 1 - testKernel / inrisk
  lambda[is.na(lambda)] <- 1
  lambda[lambda > 1 | lambda < 0] <- 1

  S <- apply(lambda, 2, cumprod)
  surv <- S[Censor, , drop = FALSE]
  timepoints <- sort(object$y[object$censor])

  return(list("surv" = surv, "timepoints" = timepoints))
}

#' Internal prediction function for regression models
#' 
#' Provides support for predicting values underneath a regression model.
#' 
#' @param object Fitted object
#' @param testx  Testing data
#' @param ...    Additional parameters, not used.
#' 
#' @return 
#' The predicted object
#' 
#' @noRd
predict_orthoDr_reg <- function(object, testx, ...) {
  # transform the covariates into the same scale
  x <- object$x
  xscale <- apply(x, 2, sd)
  xmean <- apply(x, 2, mean)
  X <- scale(x)

  testX <- sweep(testx, 2, xmean, FUN = "-")
  testX <- sweep(testX, 2, xscale, FUN = "/")

  XB <- X %*% object$B
  testXB <- testX %*% object$B

  XBscale <- apply(XB, 2, sd)
  testXB <- sweep(testXB, 2, XBscale, FUN = "/")

  XB <- scale(XB) / object$bw / sqrt(2)
  testXB <- testXB / object$bw / sqrt(2)

  # kernel matrix between x and testx
  testKernel <- KernelDist_cross(testXB, XB)

  # this method does not deal with ties. I need to fix this later on

  pred <- apply(testKernel, 1, function(w, x) weighted.mean(x, w), object$y)

  return(list("pred" = pred))
}

#' Internal prediction function for personalized dose models
#' 
#' Provides support for predicting personalized dose values.
#' 
#' @param object Fitted object
#' @param testx  Testing data
#' @param ...    Additional parameters, not used.
#' 
#' @return 
#' The predicted object
#' 
#' @noRd
predict_orthoDr_pdose <- function(object, testx, ...) {
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.matrix(testx) || !is.numeric(testx)) stop("testx must be a numercial matrix")

  if (class(object)[2] != "fit") {
    stop("This is not an orthoDr fitted object")
  }

  if (!object$keep.data) {
    stop("Need the original data for prediction. Please specify keep.data = TRUE in model fitting.")
  }

  if (class(object)[3] == "pdose") {
    pred <- dosepred(object$B, object$x, testx, object$bw, object$W)
  }

  return(list("pred" = pred))
}
