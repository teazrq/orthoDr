#' Semiparametric dimension reduction method from Ma & Zhu (2012).
#' 
#' Performs the semiparametric dimension reduction method associated with
#' Ma & Zhu (2012).
#' 
#' @param x         A `matrix` or `data.frame` for features (continous only). 
#'                  The algorithm will not scale the columns to unit variance
#' @param y         A `vector` of continuous outcome
#' @param method    Dimension reduction methods (semi-): `"sir"`, `"save"`, 
#'                  `"phd"`, `"local"` or `"seff"`. Currently only 
#'                  `"sir"` and `"phd"` are available.
#' @param ndr       The number of directions
#' @param B.initial Initial `B` values. If specified, must be a matrix with
#'                  `ncol(x)` rows and `ndr` columns. Will be processed by
#'                  Gram-Schmidt if not orthogonal. If the initial value is
#'                  not given, three initial values (`"sir"`, `"save"` and
#'                  `"phd"`) using the traditional method will be tested. The
#'                  one with smallest l2 norm of the estimating equation will
#'                  be used.
#' @param bw        A Kernel bandwidth, assuming each variables have unit
#'                  variance
#' @param keep.data Should the original data be kept for prediction. Default
#'                  is `FALSE`.
#' @param control   A list of tuning variables for optimization. 
#'                  `epsilon` is the size for numerically approximating the
#'                  gradient. For others, see Wen and Yin (2013).
#' @param maxitr    Maximum number of iterations
#' @param verbose   Should information be displayed
#' @param ncore     Number of cores for parallel computing. 
#'                  The default is the maximum number of threads.
#' @return 
#' 
#' A `orthoDr` object consisting of `list` with named elements:
#' 
#' \item{B}{The optimal `B` value}
#' \item{fn}{The final functional value}
#' \item{itr}{The number of iterations}
#' \item{converge}{convergence code}
#' 
#' @export
#' 
#' @references
#' Ma, Y., & Zhu, L. (2012). A semiparametric approach to dimension reduction. Journal of the American Statistical Association, 107(497), 168-179.
#' DOI: \doi{10.1080/01621459.2011.646925}
#' 
#' Ma, Y., & Zhu, L. (2013). Efficient estimation in sufficient dimension reduction. Annals of statistics, 41(1), 250.
#' DOI: \doi{10.1214/12-AOS1072}
#' 
#' @examples
#' # generate some regression data
#' set.seed(1)
#' N <- 100
#' P <- 4
#' dataX <- matrix(rnorm(N * P), N, P)
#' Y <- -1 + dataX[, 1] + rnorm(N)
#'
#' # fit the semi-sir model
#' orthoDr_reg(dataX, Y, ndr = 1, method = "sir")
#'
#' # fit the semi-phd model
#' Y <- -1 + dataX[, 1]^2 + rnorm(N)
#' orthoDr_reg(dataX, Y, ndr = 1, method = "phd")
orthoDr_reg <- function(x, y, method = "sir", ndr = 2,
                        B.initial = NULL, bw = NULL, keep.data = FALSE,
                        control = list(), maxitr = 500, verbose = FALSE, ncore = 0) {
  if (!is.matrix(x)) stop("x must be a matrix")
  if (!is.numeric(x)) stop("x must be numerical")
  if (nrow(x) != length(y)) stop("Number of observations do not match")

  # check tuning parameters
  control <- control.check(control)
  match.arg(method, c("sir", "save", "phd", "local", "seff"))
  # match.arg(method, c("sir", "phd"))

  # this is temporary
  # if (method %in% c("save"))
  #  stop("save and seff are currently unaviable.")

  ndr <- max(1, ndr)
  if (ndr > 4) warning("ndr > 3 is not recommended")
  ndr <- min(ndr, ncol(x))

  N <- nrow(x)
  P <- ncol(x)

  if (is.null(bw)) {
    bw <- silverman(ndr, N)
  }

  # scale y
  Y <- as.matrix(scale(y) / N^(-1 / (ndr + 5)) / sqrt(2))

  # scale x
  X <- x

  if (is.null(B.initial)) {
    B.initial <- initB(X, Y, ndr, bw, method, ncore)
  } else {
    if (!is.matrix(B.initial)) stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) | ndr != ncol(B.initial)) stop("Dimension of B.initial is not correct")
    B.initial <- gramSchmidt(B.initial)$Q
  }

  # start to fit the model
  pre <- Sys.time()

  fit <- reg_solver(
    method, B.initial, X, Y, bw,
    control$rho, control$eta, control$gamma, control$tau, control$epsilon,
    control$btol, control$ftol, control$gtol, maxitr, verbose, ncore
  )

  if (verbose > 0) {
    cat(paste("Total time: ", round(as.numeric(Sys.time() - pre, units = "secs"), digits = 2), " secs\n", sep = ""))
  }

  fit$method <- method
  fit$keep.data <- keep.data

  if (keep.data) {
    fit[["x"]] <- x
    fit[["y"]] <- y
    fit[["bw"]] <- bw
  }

  class(fit) <- c("orthoDr", "fit", "reg")

  return(fit)
}




# switch function
#' @title reg_solve
#' @name reg_solve
#' @description regression solver switch function
#' @keywords internal
#'
reg_solver <- function(method, ...) {
  if (method == "sir") {
    fit <- sir_solver(...)
  }

  if (method == "save") {
    fit <- save_solver(...)
  }

  if (method == "phd") {
    fit <- phd_solver(...)
  }

  if (method == "seff") {
    fit <- seff_solver(...)
  }

  if (method == "local") {
    fit <- local_solver(...)
  }

  return(fit)
}


# initiation function
#' @title reg_init
#' @name reg_init
#' @description regression initiation function to get better initial value
#' @keywords internal
#'
reg_init <- function(method, ...) {
  if (method == "sir") {
    f <- sir_init(...)
  }

  if (method == "save") {
    f <- save_init(...)
  }

  if (method == "phd") {
    f <- phd_init(...)
  }

  if (method == "seff") {
    f <- seff_init(...)
  }

  if (method == "local") {
    f <- local_f(...)
  }

  return(f)
}

# initiation function for B
#' @title initB
#' @name initB
#' @description regression initiation function
#' @keywords internal
#'
initB <- function(x, y, ndr, bw, method, ncore) {
  B1 <- gramSchmidt(dr(y ~ x, method = "sir")$evectors[, 1:ndr, drop = FALSE])$Q
  B2 <- gramSchmidt(dr(y ~ x, method = "save")$evectors[, 1:ndr, drop = FALSE])$Q
  B3 <- gramSchmidt(dr(y ~ x, method = "phd")$evectors[, 1:ndr, drop = FALSE])$Q

  values <- c(
    reg_init(method, B1, x, y, bw, ncore),
    reg_init(method, B2, x, y, bw, ncore),
    reg_init(method, B3, x, y, bw, ncore)
  )

  return(list(B1, B2, B3)[[which.min(values)]])
}
