#' @title predict.orthoDr
#' @name predict.orthoDr
#' @export
#' @description The prediction function for orthoDr fitted models
#' @param object A fitted orthoDr object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object. Survival predictions contain one row per
#'   observed event. Events tied at the same time have the same post-event
#'   probability and repeated time labels.
#' @details Regression and survival predictions use the fitted directions in
#'   original feature units and the training projection scales. Legacy forward
#'   survival objects without a stored bandwidth require the original bandwidth
#'   to be restored or the model to be refitted. If all Gaussian weights for a
#'   survival prediction underflow to zero, a warning is issued and that curve
#'   contains \code{NA}.
#' @examples
#' # generate some survival data
#' N = 100; P = 4; dataX = matrix(rnorm(N*P), N, P)
#' Y = exp(-1 + dataX[,1] + rnorm(N))
#' Censor = rbinom(N, 1, 0.8)
#'
#' # fit the model with keep.data = TRUE
#' orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 1,
#'                            method = "dm", keep.data = TRUE)
#'
#' #predict 10 new observations
#' predict(orthoDr.fit, matrix(rnorm(10*P), 10, P))
#'
#' # generate some personalized dose scenario
#'
#' exampleset <- function(size,ncov){
#'
#'  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
#'  A = runif(size,0,2)
#'
#'  Edr =  as.matrix(c(0.5,-0.5))
#'
#'  D_opt = X %*% Edr + 1
#'
#'  mu = 2 + 0.5*(X %*% Edr) - 7*abs(D_opt-A)
#'
#'  R = rnorm(length(mu),mu,1)
#'
#'  R = R - min(R)
#'
#'  datainfo = list(X=X,A=A,R=R,D_opt=D_opt,mu=mu)
#'  return(datainfo)
#' }
#'
#' # generate data
#'
#' set.seed(123)
#' n = 150
#' p = 2
#' ndr =1
#' train = exampleset(n,p)
#' test = exampleset(500,p)
#'
#' # the direct learning method
#' orthofit = orthoDr_pdose(train$X, train$A, train$R, ndr = ndr, lambda = 0.1,
#'                       method = "direct", K = as.integer(sqrt(n)), keep.data = TRUE,
#'                       maxitr = 150, verbose = FALSE, ncore = 2)
#'
#' predict(orthofit,test$X)
#'
#' # the pseudo direct learning method
#' orthofit = orthoDr_pdose(train$X, train$A, train$R, ndr = ndr, lambda = seq(0.1,0.2,0.01),
#'                       method = "pseudo_direct", K = as.integer(sqrt(n)), keep.data = TRUE,
#'                       maxitr = 150, verbose = FALSE, ncore = 2)
#'
#' predict(orthofit,test$X)

predict.orthoDr <- function(object, testx, ...)
{
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.matrix(testx) || !is.numeric(testx)) stop("testx must be a numerical matrix")

  if (class(object)[2] !="fit")
    stop("This is not an orthoDr fitted object")

  if (!object$keep.data)
    stop("Need the original data for prediction. Please specify keep.data = TRUE in model fitting.")

  # predict survival functions on the testing data

  if (inherits(object, "surv"))
    pred = predict_orthoDr_surv(object, testx, ...)

  # predict regression outcome on the testing data

  if (inherits(object, "reg"))
    pred = predict_orthoDr_reg(object, testx, ...)

  # predict rewards and optimal dose on the testing data

  if (class(object)[3] == "pdose")
    pred = predict_orthoDr_pdose(object, testx, ...)

  class(pred) <- c("orthoDr", "predict", class(object)[3])

  return(pred)
}

#' Internal prediction function for survival models
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#' @keywords internal
#' @noRd
predict_orthoDr_surv <- function(object, testx, ...)
{
  testKernel = prediction_kernel(object, testx)
  ord = order(object$y)
  times = object$y[ord]
  Censor = (object$censor[ord] == 1)
  testKernel = testKernel[, ord, drop = FALSE]

  # Complete risk sets include every observation at the current time.
  blocks = cumsum(c(TRUE, diff(times) != 0))
  blockWeights = rowsum(t(testKernel), blocks, reorder = FALSE)
  deaths = rowsum(t(sweep(testKernel, 2, Censor, FUN = "*")),
                  blocks, reorder = FALSE)
  S = matrix(NA_real_, nrow(blockWeights), ncol(blockWeights))
  underflow = colSums(blockWeights) == 0
  if (any(underflow))
    warning("Gaussian prediction weights underflowed to zero; returning NA survival probabilities for affected test observations.")

  for (j in which(!underflow)) {
    # Reverse sums retain representable small risk masses late in follow-up.
    risk = rev(cumsum(rev(blockWeights[, j])))
    d = deaths[, j]
    if (any(!is.finite(risk) | !is.finite(d) | risk < 0 | d < 0))
      stop("Invalid weighted survival risk or death mass")
    excess = d - risk
    # Independently summed masses may differ by at most 64 eps relatively.
    tolerance = 64 * .Machine$double.eps * pmax(abs(d), abs(risk))
    if (any(excess > tolerance))
      stop("Weighted survival death mass exceeds its risk mass")
    d = pmin(d, risk)
    factors = rep(1, length(risk))
    positive = d > 0
    factors[positive] = 1 - d[positive] / risk[positive]
    S[, j] = cumprod(factors)
  }

  surv = S[blocks[Censor], , drop = FALSE]
  timepoints = times[Censor]

  return(list("surv" = surv, "timepoints" = timepoints))
}

#' Internal prediction function for regression models
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#' @keywords internal
#' @noRd
predict_orthoDr_reg <- function(object, testx, ...)
{
  testKernel = prediction_kernel(object, testx)

  pred = apply(testKernel, 1, function(w, x) weighted.mean(x, w), object$y)

  return(list("pred" = pred))
}

#' Gaussian prediction kernel in the returned raw-feature basis
#' @keywords internal
#' @noRd
prediction_kernel <- function(object, testx)
{
  if (is.null(object$bw))
    stop("This fitted object has no bandwidth. Supply the original bandwidth as object$bw or refit the model.")
  if (length(object$bw) != 1L || !is.finite(object$bw) || object$bw <= 0)
    stop("The fitted bandwidth must be a positive finite scalar")

  XB = object$x %*% object$B
  testXB = testx %*% object$B
  XBmean = colMeans(XB)
  XBscale = apply(XB, 2, sd)
  if (any(!is.finite(XBscale) | XBscale <= 0))
    stop("Training projection scales must be positive and finite for prediction")

  denominator = XBscale * object$bw * sqrt(2)
  XB = sweep(sweep(XB, 2, XBmean, FUN = "-"), 2, denominator, FUN = "/")
  testXB = sweep(sweep(testXB, 2, XBmean, FUN = "-"), 2, denominator, FUN = "/")
  weights = KernelDist_cross(testXB, XB)
  if (any(!is.finite(weights) | weights < 0))
    stop("Prediction kernel weights must be finite and nonnegative")
  weights
}

#' Internal prediction function for personalized dose models
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#' @keywords internal
#' @noRd
predict_orthoDr_pdose <- function(object, testx, ...)
{
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.matrix(testx) || !is.numeric(testx)) stop("testx must be a numerical matrix")

  if (class(object)[2] !="fit")
    stop("This is not an orthoDr fitted object")

  if (!object$keep.data)
    stop("Need the original data for prediction. Please specify keep.data = TRUE in model fitting.")

  if (class(object)[3] == "pdose")
      pred = dosepred(object$B, object$x, testx, object$bw, object$W)

  return(list("pred" = pred))
}


