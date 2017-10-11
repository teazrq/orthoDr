#' @title IR-CP model
#' @name orthoDr_surv
#' @description The counting process based semiparametric dimention reduction (IR-CP) model for right censored survival outcome.
#' @param x A matrix or data.frame for features (continous only). The algorithm will not scale the columns to unit variance
#' @param y A vector of observed time
#' @param censor A vector of censoring indicator
#' @param ndr The number of directions
#' @param B.initial Initial \code{B} values. Will use the counting process based SIR model \link[orthoDr]{CP_SIR} as the initial if leaving as \code{NULL}.
#' If specified, must be a matrix with \code{ncol(x)} rows and \code{ndr} columns. Will be processed by Gram-Schmidt if not orthogonal
#' @param bw A Kernel bandwith, assuming each variables have unit variance
#' @param control A list of tuning varaibles for optimization. \code{epsilon} is the size for numerically appriximating the gradient. For others, see Wen and Yin (2013).
#' @param maxitr Maximum number of iterations
#' @param verbose Should information be displayed
#' @return A \code{orthoDr} object; a list consisting of
#' \item{B}{The optimal \code{B} value}
#' \item{fn}{The final funtional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' @references Sun, Q., Zhu, R., Wang T. and Zeng D. "Counting Process Based Dimension Reduction Method for Censored Outcomes." (2017)
#' \url{https://arxiv.org/abs/1704.05046} .
#' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434.
#' DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
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
#'
#' # fit the model
#' orthoDr.fit = orthoDr_surv(dataX, Y, Censor, ndr = 2)
#'
#' # compare with the true direction
#' c(distance(failEDR, orthoDr.fit$B, "dist"),
#'   distance(failEDR, orthoDr.fit$B, "trace"),
#'   distance(failEDR, orthoDr.fit$B, "canonical", dataX))


orthoDr_surv <- function(x, y, censor, ndr = 2, B.initial = NULL, bw = NULL,
                         control = list(), maxitr = 500, verbose = FALSE)
{
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y) | nrow(x) != length(censor)) stop("Number of observations do not match")

  if (is.null(bw))
    bw = silverman(ndr, nrow(x))

  if (is.null(B.initial))
  {
    B.initial = CP_SIR(x, y, censor)$vectors[,1:ndr, drop = FALSE]
  }else{
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) | ndr != ncol(B.initial)) stop("Dimention of B.initial is not correct")
  }


  # check tuning parameters

  if (!is.list(control))
    stop("control must be a list of tuning parameters")

  if (is.null(control$rho))
  {
    control$rho = 1e-4
  }else if(control$rho < 0) control$rho = 1e-4

  if (is.null(control$eta))
  {
    control$eta = 0.2
  }else if(control$eta < 0) control$eta = 0.2

  if (is.null(control$gamma))
  {
    control$gamma = 0.85
  }else if(control$gamma < 0) control$gamma = 0.85

  if (is.null(control$tau))
  {
    control$tau = 1e-3
  }else if(control$tau < 0) control$tau = 1e-3

  if (is.null(control$epsilon))
  {
    control$epsilon = 1e-6
  }else if(control$epsilon < 0) control$epsilon = 1e-6

  if (is.null(control$btol))
  {
    control$btol = 1e-6
  }else if(control$btol < 0) control$btol = 1e-6

  if (is.null(control$ftol))
  {
    control$ftol = 1e-6
  }else if(control$ftol < 0) control$ftol = 1e-6

  if (is.null(control$gtol))
  {
    control$gtol = 1e-6
  }else if(control$gtol < 0) control$gtol = 1e-6

	# center matrix X, but do not scale

  B.initial = gramSchmidt(B.initial)$Q

  N = nrow(x)
  P = ncol(x)

  Yorder = order(y)
  X = scale(x, scale = FALSE)
  X = X / sd(as.vector(X))
  X = X[Yorder, ]
  Y = seq(0, 1, length.out = length(y))
  C = censor[Yorder]
  Fail.Ind = which(C==1)

  # calculate some useful stuff

  timepoints = Y[C == 1]
  nFail = length(Fail.Ind)

  # inRisk set at all failure times

  inRisk = matrix(NA, N, nFail)

  for (j in 1:nFail)
    inRisk[, j] = as.matrix(Y >= timepoints[j])

  # E[X | dN(t) = 1, Y(t) = 1] - E[X | dN(t) = 0, Y(t) = 1] at all failure times

  kernel.y = exp(-(as.matrix(dist(Y, method = "euclidean"))/ silverman(1, length(y)))^2)

  Phit = matrix(0, P, nFail)

  for (j in 1:nFail)
    Phit[,j] = as.matrix(apply(X, 2, weighted.mean, w = C*kernel.y[, Fail.Ind[j]]) - colMeans(X[Fail.Ind[j]:N, ,drop = FALSE]))

  # start to fit the model

  pre = Sys.time()
  fit = surv_solver(B.initial, X, Phit, inRisk, bw, Fail.Ind,
                    control$rho, control$eta, control$gamma, control$tau, control$epsilon,
                    control$btol, control$ftol, control$gtol, maxitr, verbose)
  if (verbose > 0)
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))

	class(fit) <- "orthoDr"

	return(fit)
}
