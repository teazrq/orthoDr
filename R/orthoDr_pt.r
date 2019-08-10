#' @title orthoDr_pm model
#' @name orthoDr_pm
#' @description The "Direct Learning & Pseudo-direct Learning" Method for personalized medicine.
#' @param x A matrix or data.frame for features (continuous only).
#' @param a A vector of observed dose
#' @param r A vector of observed reward
#' @param ndr A dimension structure
#' @param B.initial Initial \code{B} values. Will use the counting process based SIR model \link[orthoDr]{CP_SIR} as the initial if leaving as \code{NULL}.
#' If specified, must be a matrix with \code{ncol(x)} rows and \code{ndr} columns. Will be processed by Gram-Schmidt if not orthogonal
#' @param bw A Kernel bandwidth, assuming each variables have unit variance
#' @param lambda A GCV penalty for the kernel ridge regression
#' @param K A number of grids in the range of dose
#' @param method A method the user will implement
#' @param keep.data Should the original data be kept for prediction
#' @param control A list of tuning variables for optimization. \code{epsilon} is the size for numerically appriximating the gradient. For others, see Wen and Yin (2013).
#' @param maxitr Maximum number of iterations
#' @param ncore the number of cores for parallel computing
#' @param verbose Should information be displayed
#' @return A \code{orthoDr} object; a list consisting of
#' \item{B}{The optimal \code{B} value}
#' \item{fn}{The final functional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' @references Zhou, W., Zhu, R. "A Parsimonious Personalized Dose Model vis Dimension Reduction." (2018+)
#' \url{https://arxiv.org/abs/1802.06156}.
#' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434.
#' DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
#' @examples
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
#' dose = predict(orthofit,test$X)
#'
#` # compare with the optimal dose
#' dosedistance = mean((test$D_opt-dose$pred)^2)
#' print(dosedistance)
#' 
#' the pseudo direct learning method 					  
#' orthofit = orthoDr_pm(train$X,train$A,train$R,ndr = ndr,lambda = seq(0.1,0.2,0.01),
#'                       method = "pseudo_direct", keep.data = T)
#'                     
#'dose = predict(orthofit,test$X)
#'
#' compare with the optimal dose
#'dosedistance = mean((test$D_opt-dose$pred)^2)
#'print(dosedistance)

orthoDr_pm <- function(x, a, r, ndr = ndr, B.initial = NULL, bw = NULL, lambda = 0.1,
                       K = sqrt(length(r)), method = c("direct","pseudo_direct"),
                       keep.data = FALSE, control = list(), maxitr = 500, verbose = FALSE, ncore = 0)
{
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(r) | nrow(x) != length(a)) stop("Number of observations do not match")

  if (is.null(bw))
    bw = silverman(ndr, nrow(x))
  if (is.null(B.initial))
  {
    n= nrow(x)
    p = ncol(x)
    B.initial = pSave(x, a, r, ndr = ndr)
  }else{
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) | ndr != ncol(B.initial)) stop("Dimention of B.initial is not correct")
  }

  # check tuning parameters
  control = control.check(control)

  B.initial = gramSchmidt(B.initial)$Q

  N = nrow(x)
  P = ncol(x)
  X = x
  
  # standerdize
  a_center = mean(a)
  a_scale = sd(a)
  a_scale_bw = a/a_scale/bw

  
  cdose= seq(min(a), max(a), length.out = K)
  
  cdose_scale = cdose/sd(cdose)/bw
  
  A.dist <- matrix(NA, nrow(x), K)
  for (k in 1:K)
  {
    A.dist[, k] = exp(-((a_scale_bw - cdose_scale[k]))^2)
  }
  
  if (method == "direct")
  {
    
    pre = Sys.time()
    fit = direct_pt_solver(B.initial, X, a, A.dist, cdose, r, lambda, bw,
                           control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                           control$btol, control$ftol, control$gtol, maxitr, verbose, ncore)
    if (verbose > 0)
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
  }

  if(method == "pseudo_direct")
  {

    pre = Sys.time()
    fit = semi_pt_solver(B.initial, X, r, a, A.dist, cdose, lambda, bw,
                         control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                         control$btol, control$ftol, control$gtol, maxitr, verbose,ncore)
    if (verbose > 0)
      cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))

  }

  fit$method = method
  fit$keep.data = keep.data

  if (keep.data)
  {
    fit[['x']] = x
    fit[['a']] = a
    fit[['r']] = r
    fit[['bw']] = bw
  }

  class(fit) <- c("orthoDr", "fit", "personalized_treatment", method)

  return(fit)
}

