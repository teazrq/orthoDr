#' @title predict.orthoDr
#' @name predict.orthoDr
#' @description The prediction function for orthoDr fitted models
#' @param object A fitted orthoDr object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#' @examples
#' # generate some survival data
#' N = 100; P = 4; dataX = matrix(rnorm(N*P), N, P)
#' Y = exp(-1 + dataX[,1] + rnorm(N))
#' Censor = rbinom(N, 1, 0.8)
#'
#' # fit the model with keep.data = TRUE
#' orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 1, method = "dm", keep.data = TRUE)
#'
#' #predict 10 new observations
#' predict(orthoDr.fit, matrix(rnorm(10*P), 10, P))
#'
#' # generate some personalized dose scenario
#' 
#' Scenario <- function(size,ncov)
#' {
#'  set.seed(as.integer((as.double(Sys.time())*100+Sys.getpid()) %% 2^14) )
#'  X = matrix(runif(size*ncov,-1,1),ncol=ncov)
#'  A = runif(size,0,2)
#'
#'  Edr = as.matrix(cbind(c(1, 0.5,0, 0, -0.5, 0, 0, 0,rep(0,2)),
#'                        c(0.5, 0, 0.5, -0.5, 1,0,0,0,rep(0, 2))))
#'
#'  D_opt = sin(X %*% Edr[,2] * X %*% Edr[,1]) + (3/4)*(X %*% Edr[,1])/ (5 + (X %*% Edr[,2] + 4)^2) + 1
#'
#'  mu = 7 + 0.5*(X %*% Edr[,1])^2 + 1*X %*% Edr[,2] - 13*abs(D_opt-A)
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
#' n = 400
#' p = 10
#' ndr =2
#' train = Scenario1(n,p)
#' test = Scenario1(500,p)
#'
#' # the pseudo direct learning method
#'  orthoDr_pm(train$X,train$A,train$R,ndr =ndr,lambda = seq(0.1,0.2,0.01), 
#'             method = "direct",keep.data = T)
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

  if (class(object)[3] == "survival")
    pred = predict_orthoDr_surv(object, testx, ...)

  # predict regression outcome on the testing data

  if (class(object)[3] == "regression")
    pred = predict_orthoDr_reg(object, testx, ...)

  # predict rewards and optimal dose on the testing data

  if (class(object)[3] == "personalized_treatment")
    pred = predict_orthoDr_pt(object, testx, ...)

  class(pred) <- c("orthoDr", "predict", class(object)[3])

  return(pred)
}


#' @title predict_orthoDr_surv
#' @name predict_orthoDr_surv
#' @description Internal prediction function for survival models
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#'
predict_orthoDr_surv <- function(object, testx, ...)
{
  # transform the covariates into the same scale
  x = object$x
  xscale = apply(x, 2, sd)
  xmean = apply(x, 2, mean)
  X = scale(x)

  testX = sweep(testx, 2, xmean, FUN = "-")
  testX = sweep(testX, 2, xscale, FUN = "/")

  XB = X %*% object$B
  testXB = testX %*% object$B

  XBscale = apply(XB, 2, sd)
  testXB = sweep(testXB, 2, XBscale, FUN = "/")

  XB = scale(XB)/object$bw/sqrt(2)
  testXB = testXB/object$bw/sqrt(2)

  # kernel matrix between x and testx
  testKernel = KernelDist_cross(testXB, XB)

  # this method does not deal with ties. I need to fix this later on

  testKernel = testKernel[, order(object$y), drop = FALSE]
  Censor = (object$censor[order(object$y)] == 1)

  inrisk = apply(testKernel, 1, cumsum)
  totalweights = inrisk[nrow(X), , drop = FALSE]
  inrisk = sweep(-inrisk, 2, totalweights, FUN = "+")

  testKernel = sweep(testKernel, 2, Censor, FUN = "*" )
  testKernel = t(testKernel)

  lambda = 1 - testKernel / inrisk
  lambda[is.na(lambda)] = 1
  lambda[lambda > 1 | lambda < 0] = 1

  S = apply(lambda, 2, cumprod)
  surv = S[Censor, , drop = FALSE]
  timepoints = sort(object$y[object$censor])

  return(list("surv" = surv, "timepoints" = timepoints))
}


#' @title predict_orthoDr_reg
#' @name predict_orthoDr_reg
#' @description Internal prediction function for regression models
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object
#'
predict_orthoDr_reg <- function(object, testx, ...)
{
  # transform the covariates into the same scale
  x = object$x
  xscale = apply(x, 2, sd)
  xmean = apply(x, 2, mean)
  X = scale(x)

  testX = sweep(testx, 2, xmean, FUN = "-")
  testX = sweep(testX, 2, xscale, FUN = "/")

  XB = X %*% object$B
  testXB = testX %*% object$B

  XBscale = apply(XB, 2, sd)
  testXB = sweep(testXB, 2, XBscale, FUN = "/")

  XB = scale(XB)/object$bw/sqrt(2)
  testXB = testXB/object$bw/sqrt(2)

  # kernel matrix between x and testx
  testKernel = KernelDist_cross(testXB, XB)

  # this method does not deal with ties. I need to fix this later on

  pred = apply(testKernel, 1, function(w, x) weighted.mean(x, w), object$y)

  return(list("pred" = pred))
}


#' @title predict_orthoDr_pt
#' @name predict_orthoDr_pt
#' @description Internal prediction function for personalized treatment models
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object


predict_orthoDr_pt <- function(object, testx, ...)
{
  # check test data
  if (missing(testx)) stop("testx is missing")
  if (!is.matrix(testx) || !is.numeric(testx)) stop("testx must be a numercial matrix")

  if (class(object)[2] !="fit")
    stop("This is not an orthoDr fitted object")

  if (!object$keep.data)
    stop("Need the original data for prediction. Please specify keep.data = TRUE in model fitting.")


  if (class(object)[3] == "personalized_treatment")
    if  (class(object)[4] == "direct"){
      pred = predict_orthoDr_pt_direct(object, testx, ...)
    }
  if  (class(object)[4] == "pseudo_direct"){
      pred = predict_orthoDr_pt_pseudo_direct(object, testx, ...)
  }

  class(pred) <- c("orthoDr", "predict", class(object)[3],class(object)[4])

  return(pred)
}

#' @title predict_orthoDr_pt_direct
#' @name predict_orthoDr_pt_direct_kernel
#' @description Internal prediction function for direct personalized dose model
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object

predict_orthoDr_pt_direct<- function(object, testx, ...)
{

  pred = Dosepred(object$B, object$x, testx, object$bw, object$W)
  return(list("pred" = pred))

}

#' @title predict_orthoDr_pt_pseudo_direct
#' @name predict_orthoDr_pt_semi_svm
#' @description Internal prediction function for pseudo_direct personalized dose model
#' @keywords internal
#' @param object fitted object
#' @param testx Testing data
#' @param ... ...
#' @return The predicted object

predict_orthoDr_pt_pseudo_direct<- function(object, testx, ...)
{

  v = list()
  v$X = object$x
  v$A = object$a
  v$R = object$r

  index = which(v$R > quantile(v$R,0.65))
  model_nopen  = svm(x = (v$X %*% as.matrix(object$B))[index,], y = v$A[index], w= v$R[index], type="eps-regression",
                     epsilon = 0.15, scale=F)

  pred = predict(model_nopen, testx %*% as.matrix(object$B))
  return(list("pred" = pred))

}


