#' Semiparametric Dimension Reduction for Regression
#'
#' @description
#' Fits a semiparametric sufficient dimension reduction model for continuous
#' outcomes using the methods of Ma & Zhu (2012, 2013). It estimates a
#' direction matrix \code{B} subject to \eqn{B^\top B = I}; \eqn{xB} gives
#' the reduced predictors.
#'
#' Kernel calculations divide each coordinate of \eqn{B^\top x} by its sample
#' standard deviation and by \code{bw * sqrt(2)}. SIR, PHD, LOCAL, and SEFF
#' center the predictors and divide by one pooled sample standard deviation
#' for optimization. A common change of units preserves these scaled objectives.
#' SAVE uses its supplied-unit estimating equation and divides its squared norm
#' by \eqn{\max(s^4, s^8)}, where \eqn{s} is the pooled sample standard deviation.
#' A common change of units can change the SAVE fit because the equation
#' includes the identity minus a conditional covariance matrix. Rescaling
#' individual predictors can affect all five methods.
#' With a very narrow kernel, a finite-difference step can move a projected
#' observation by several bandwidths. If using such a bandwidth with SIR or
#' SAVE, check whether a smaller \code{control$epsilon} changes the fitted
#' directions or held-out predictions.
#' With \code{keep.data = TRUE}, the original predictors
#' are retained, and prediction accepts new data in the supplied coordinates.
#' If predictors are standardized before fitting, apply the same training
#' transformation to prediction data. The outcome \code{y} is centered and
#' scaled internally for estimation; mean predictions use the original
#' outcome units.
#'
#' @name orthoDr_reg
#' @export
#'
#' @param x A numeric matrix of predictors (n observations \eqn{\times} p
#'   features). A single predictor must be finite and varying.
#' @param y A numeric vector or one-column matrix of finite, varying continuous
#'   outcomes, one per observation.
#' @param method The dimension reduction method. One of \code{"sir"},
#'   \code{"save"}, \code{"phd"}, \code{"local"}, or \code{"seff"}.
#'   See Details for method-specific behaviour.
#' @param ndr Number of directions to estimate (default 2). The fit uses
#'   at least one direction and at most \code{ncol(x)}. Values > 4 produce a warning.
#' @param B.initial Initial \eqn{B} matrix (\eqn{p \times ndr}), or \code{NULL}
#'   (default). With one predictor, the default is \code{matrix(1)}.
#'   With multiple predictors, the default selects the SIR, SAVE or PHD
#'   starting matrix with the smallest squared estimating-equation norm.
#'   For \code{ndr = 2}, it also considers the leading
#'   SIR and SAVE directions together when they are linearly independent.
#'   For \code{method = "local"} with \code{ndr = 1}, it instead fits from
#'   each of the three starting matrices and selects the result with the
#'   largest nonparametric \eqn{R^2} for the Nadaraya-Watson regression of
#'   \eqn{Y} on \eqn{B^T X}.
#' @param bw X-kernel bandwidth. \code{NULL} (default) uses Silverman's rule
#'   of thumb. Internal Y-kernel bandwidths are derived from the sample size
#'   and structural dimension; they are not user-tunable.
#' @param keep.data Logical; if \code{TRUE}, the original data is stored in the
#'   returned object for use by \code{\link{predict}}. Default \code{FALSE}.
#' @param control A named list of optimizer tuning parameters. Recognised
#'   entries and their defaults:
#'   \describe{
#'     \item{\code{rho}}{Coefficient in the line-search sufficient-decrease
#'       condition. Default: \code{1e-4}.}
#'     \item{\code{eta}}{Factor for decreasing the step size in backtracking
#'       line search. Default: \code{0.2}.}
#'     \item{\code{gamma}}{Weight for updating the nonmonotone line-search
#'       reference value. Default: \code{0.85}.}
#'     \item{\code{tau}}{Initial step size. Default: \code{1e-3}.}
#'     \item{\code{epsilon}}{Finite-difference perturbation.
#'       Default: \code{1e-6}.}
#'     \item{\code{btol}}{Parameter-change tolerance. Default: \code{1e-6}.}
#'     \item{\code{ftol}}{Objective-change tolerance. Default: \code{1e-6}.}
#'     \item{\code{gtol}}{Gradient-norm tolerance. Default: \code{1e-6}.}
#'   }
#' @param maxitr Maximum number of iterations (default 500). Clamped to
#'   \eqn{\ge} 5 internally.
#' @param verbose Integer; set > 0 to print iteration progress and total
#'   elapsed time. Default 0.
#' @param ncore Integer; number of CPU threads for parallel gradient
#'   approximation. 0 (default) uses all available threads.
#'
#' @details
#' \subsection{Methods}{
#'   \describe{
#'     \item{SIR}{Semiparametric sliced inverse regression (Ma & Zhu 2012).}
#'     \item{SAVE}{Semiparametric sliced average variance estimation.}
#'     \item{PHD}{Semiparametric principal Hessian directions. The estimating
#'       equation uses outer products of the centered predictors.}
#'     \item{Local}{Locally efficient estimator (Ma & Zhu 2013,
#'       Section~3.1) with a normal posited model for the conditional
#'       density. This yields a mean-regression estimator, which captures
#'       dependence of \eqn{Y} on \eqn{\beta^\top x} through the conditional
#'       mean only, not higher moments. For \eqn{d = 1 < p} with default
#'       initialisation, a multi-start strategy runs LOCAL from SIR, SAVE,
#'       and PHD starters and selects the result that maximises the
#'       nonparametric \eqn{R^2} of the Nadaraya–Watson regression of
#'       \eqn{Y} on \eqn{\beta^\top x}.}
#'     \item{SEFF}{Semiparametric efficient estimator (Ma & Zhu 2013,
#'       Section 3.2). During optimization, it
#'       uses the initial projected covariates as a fixed kernel reference
#'       and re-evaluates the nuisance components at each candidate projection.}
#'   }
#' }
#'
#' @return An object of class \code{c("orthoDr", "fit", "reg")}, a list with
#'   components:
#'   \describe{
#'     \item{\code{B}}{Estimated \eqn{p \times ndr} direction matrix with
#'       orthonormal columns.}
#'     \item{\code{fn}}{Final objective value in the supplied predictor units.
#'       Optimizer output and stopping tolerances use the internal scaled
#'       objective. At extreme predictor scales the reported value can round
#'       to zero or overflow even when the fit converges.}
#'     \item{\code{itr}}{Number of iterations performed.}
#'     \item{\code{converge}}{Logical. \code{TRUE} indicates a small constrained
#'       gradient or small changes in both the basis and the objective, measured
#'       in the latest step or averaged over the last five steps.}
#'     \item{\code{method}}{Method used (e.g., \code{"sir"}).}
#'     \item{\code{keep.data}}{Whether original data was stored.}
#'   }
#'   If \code{keep.data = TRUE}, the list also contains \code{x}, \code{y},
#'   and \code{bw} for prediction.
#'
#' @references
#' Ma, Y. and Zhu, L. (2012). A semiparametric approach to dimension reduction.
#' \emph{Journal of the American Statistical Association}, 107(497), 168–179.
#' \doi{10.1080/01621459.2011.646925}.
#'
#' Ma, Y. and Zhu, L. (2013). Efficient estimation in sufficient dimension
#' reduction. \emph{Annals of Statistics}, 41(1), 250–268.
#' \doi{10.1214/12-AOS1072}.
#'
#' @seealso \code{\link{predict.orthoDr}}, \code{\link{ortho_optim}}
#'
#' @examples
#' set.seed(1)
#' N <- 100; P <- 4
#' X <- matrix(rnorm(N * P), N, P)
#'
#' # Linear mean model (SIR)
#' Y <- -1 + X[, 1] + rnorm(N)
#' fit_sir <- orthoDr_reg(X, Y, ndr = 1, method = "sir")
#' fit_sir$B
#'
#' # Quadratic mean model (PHD)
#' Y <- -1 + X[, 1]^2 + rnorm(N)
#' fit_phd <- orthoDr_reg(X, Y, ndr = 1, method = "phd")
#' fit_phd$B
#'
#' # Efficient estimator (SEFF)
#' fit_seff  <- orthoDr_reg(X, Y, ndr = 1, method = "seff")
#' fit_seff$B
#'
#' # Custom initial + prediction
#' B0 <- matrix(c(1, 0, 0, 0), 4, 1)
#' fit <- orthoDr_reg(X, Y, ndr = 1, method = "phd",
#'                    B.initial = B0, keep.data = TRUE)
#' predict(fit, X[1:5, ])
#'
orthoDr_reg <- function(x, y,
                        method    = c("sir", "save", "phd", "local", "seff"),
                        ndr       = 2,
                        B.initial = NULL,
                        bw        = NULL,
                        keep.data = FALSE,
                        control   = list(),
                        maxitr    = 500,
                        verbose   = 0,
                        ncore     = 0)
{
  # --- input validation ---
  if (!is.matrix(x))
    stop("x must be a matrix")
  if (!is.numeric(x))
    stop("x must be numerical")
  if (nrow(x) != length(y))
    stop("Number of observations do not match")

  method   <- match.arg(method)
  control  <- control.check(control)

  maxitr <- max(5L, as.integer(maxitr))

  ndr <- max(1L, as.integer(ndr))
  if (ndr > 4L)
    warning("ndr > 4 is not recommended")
  ndr <- min(ndr, ncol(x))

  N <- nrow(x)
  P <- ncol(x)

  if (P == 1L) {
    if (N < 2L)
      stop("One-predictor regression requires at least two observations")
    if (any(!is.finite(x)) || all(x == x[1L]))
      stop("One-predictor regression requires a finite, varying predictor")
    if (is.null(B.initial))
      B.initial <- matrix(1, 1L, 1L)
  }
  if (!(is.numeric(y) || (P > 1L && is.logical(y))) ||
      any(!is.finite(y)) || all(y == y[1L]))
    stop("Regression requires a finite, varying response")
  if (!is.null(dim(y)) && !(is.matrix(y) && ncol(y) == 1L))
    stop("The response must have one value per observation")

  if (is.null(bw))
    bw <- silverman(ndr, N)

  # --- data preparation ---
  # y is internally scaled per Ma & Zhu convention
  Y <- as.matrix(scale(y) / N^(-1 / (ndr + 5)) / sqrt(2))

  working <- reg_working_predictors(x, method)
  X <- working$x
  pooled_ratio <- 0
  objective_power <- if (method %in% c("local", "seff")) 2L else 4L
  save_weights <- c(identity_weight = 1, covariance_weight = 1, moment_scale = 1)
  if (method == "save") {
    local <- sir_local_predictors(x, working$scale)
    X <- local$x
    save <- save_objective_scale(working$scale, local$pooled_ratio)
    save_weights <- save$weights
    objective_power <- save$power
  }

  # --- initial value ---
  use_default_init <- is.null(B.initial)
  if (use_default_init) {
    candidates <- reg_candidates(working$x, Y, ndr)
  } else {
    if (!is.matrix(B.initial))
      stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) || ndr != ncol(B.initial))
      stop("Dimension of B.initial is not correct")
    B.initial <- orthonormalize_start(B.initial)
    candidates <- list(B.initial)
  }
  if (method == "sir" && sir_needs_local_geometry(X, candidates, bw)) {
    local <- sir_local_predictors(x, working$scale)
    X <- local$x
    pooled_ratio <- local$pooled_ratio
  }
  if (use_default_init && !(method == "local" && ndr == 1L))
    B.initial <- select_reg_candidate(candidates, X, Y, bw, method, ncore,
                                      pooled_ratio, save_weights)

  # --- multi-start for local method (ndr = 1 only) ---
  # The semi-parametric estimating-equation surface has multiple stationary
  # points.  The squared-norm objective is degenerate for directions where
  # B^T X is weakly related to Y: the kernel derivative collapses and the
  # estimating-function norm drops below the true root.  The objective alone
  # cannot distinguish a good solution from a spurious one.
  #
  # Strategy (ndr = 1): run LOCAL from three structured inits (SIR, SAVE,
  # PHD) and select the result that maximises the nonparametric R^2 of the
  # Nadaraya–Watson regression of Y on B^T X.  This R^2 is the same NW
  # estimator that defines the LOCAL method itself, so the criterion is
  # internally consistent and free of parametric assumptions.
  if (method == "local" && ndr == 1L && use_default_init) {
    t_start <- Sys.time()

    fit_sir <- reg_solver(
      "local", candidates[[1L]], X, Y, bw,
      control$rho, control$eta, control$gamma,
      control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol,
      maxitr, 0L, ncore
    )
    fit_save <- reg_solver(
      "local", candidates[[2L]], X, Y, bw,
      control$rho, control$eta, control$gamma,
      control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol,
      maxitr, 0L, ncore
    )
    fit_phd <- reg_solver(
      "local", candidates[[3L]], X, Y, bw,
      control$rho, control$eta, control$gamma,
      control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol,
      maxitr, 0L, ncore
    )

    # Nonparametric R^2 via the same NW estimator that LOCAL uses internally.
    # The correct direction should capture more variation in Y than a
    # spurious one, where B^T X is nearly independent of Y.
    nw_r2 <- function(B) {
      idx     <- X %*% B
      idx_scl <- idx / (stats::sd(as.vector(idx)) * bw * sqrt(2))
      K       <- kernel_weight(idx_scl, idx_scl)
      Kx      <- colSums(K)
      Ey      <- as.vector(K %*% Y / Kx)
      1 - sum((Y - Ey)^2) / sum((Y - mean(Y))^2)
    }

    r2_sir  <- nw_r2(fit_sir$B)
    r2_save <- nw_r2(fit_save$B)
    r2_phd  <- nw_r2(fit_phd$B)

    fit <- switch(which.max(c(r2_sir, r2_save, r2_phd)),
      fit_sir, fit_save, fit_phd
    )
    fit$fn <- reg_original_objective(fit$fn, working$scale, objective_power)

    fit$method    <- method
    fit$keep.data <- keep.data
    if (keep.data) { fit$x <- x; fit$y <- y; fit$bw <- bw }
    class(fit) <- c("orthoDr", "fit", "reg")

    if (verbose > 0)
      cat(sprintf("orthoDr_reg (local): multi-start (SIR/SAVE/PHD), NW-R^2 = %.3f / %.3f / %.3f, %.1f secs\n",
                  r2_sir, r2_save, r2_phd,
                  as.numeric(Sys.time() - t_start, units = "secs")))

    return(fit)
  }

  # --- fit ---
  t0 <- Sys.time()
  if (verbose > 0)
    cat(sprintf("%s optimizer diagnostics use the predictor-scale-adjusted objective.\n",
                toupper(method)))

  fit <- reg_solver(
    method, B.initial, X, Y, bw,
    control$rho, control$eta, control$gamma,
    control$tau, control$epsilon,
    control$btol, control$ftol, control$gtol,
    maxitr, verbose, ncore, pooled_ratio = pooled_ratio, save_weights = save_weights
  )
  fit$fn <- reg_original_objective(fit$fn, working$scale, objective_power)

  if (verbose > 0) {
    elapsed <- as.numeric(Sys.time() - t0, units = "secs")
    cat(sprintf("orthoDr_reg (%s): %d iterations, %.1f secs\n",
                method, fit$itr, elapsed))
  }

  # --- assemble result ---
  fit$method    <- method
  fit$keep.data <- keep.data

  if (keep.data) {
    fit$x  <- x
    fit$y  <- y
    fit$bw <- bw
  }

  class(fit) <- c("orthoDr", "fit", "reg")
  fit
}


# ---------------------------------------------------------------------------
#  Internal helpers
# ---------------------------------------------------------------------------

reg_working_predictors <- function(x, method) {
  if (any(!is.finite(x)) || nrow(x) < 2L)
    stop(toupper(method), " requires finite predictors and at least two observations")
  differences <- sweep(x, 2L, as.double(x[1L, ]), "-")
  shrink <- 1
  if (any(!is.finite(differences))) {
    # Halving first allows differences between opposite finite extremes.
    shrink <- 2
    differences <- sweep(x / shrink, 2L, x[1L, ] / shrink, "-")
  }
  peak <- max(abs(differences))
  if (peak == 0)
    stop(toupper(method), " requires predictor variation")
  unit <- differences / peak
  centered <- sweep(unit, 2L, colMeans(unit), "-")
  rms <- sqrt(sum(centered^2) / ((nrow(x) - 1) * ncol(x)))
  list(x = centered / rms, scale = c(shrink, peak, rms))
}

sir_local_predictors <- function(x, scale) {
  # A very narrow bandwidth can resolve differences erased by dividing
  # already centered coordinates. Delete only offsets that can be subtracted
  # exactly and use a binary scale before forming local differences.
  lower <- apply(x, 2L, min)
  upper <- apply(x, 2L, max)
  origin <- ifelse(lower > 0 & upper / 2 <= lower, lower,
                   ifelse(upper < 0 & lower / 2 >= upper, upper, 0))
  exponent <- max(-1074, min(1023, sum(floor(log2(scale)))))
  binary_scale <- 2^exponent
  list(x = sweep(x, 2L, origin, "-") / binary_scale,
       pooled_ratio = 1 / ((scale[2] / binary_scale) * scale[1] * scale[3]))
}

sir_needs_local_geometry <- function(x, candidates, bw) {
  epsilon <- .Machine$double.eps
  if (bw < sqrt(epsilon))
    return(TRUE)
  peak <- apply(abs(x), 2L, max)
  any(vapply(candidates, function(B) {
    projected <- x %*% B
    width <- apply(projected, 2L, stats::sd) * bw
    any((epsilon * peak) %*% abs(B) > sqrt(epsilon) * width)
  }, logical(1)))
}

save_objective_scale <- function(scale, pooled_ratio) {
  # A SAVE moment is s^2 A - s^4 D in terms of predictors divided by
  # their pooled SD s. Keep this ratio of terms while bounding both weights.
  large <- scale[2L] >= 1 / (scale[1L] * scale[3L])
  weights <- if (large) {
    c((1 / scale[2L] / scale[1L] / scale[3L])^2, 1)
  } else {
    c(1, (scale[2L] * (scale[1L] * scale[3L]))^2)
  }
  weights[2L] <- weights[2L] * pooled_ratio^2
  list(weights = setNames(c(weights, pooled_ratio^4),
                         c("identity_weight", "covariance_weight", "moment_scale")),
       power = if (large) 8L else 4L)
}

reg_original_objective <- function(fn, scale, power) {
  if (fn == 0 || !is.finite(fn))
    return(fn)
  # Combine bounded mantissas and binary exponents. Neither s^power nor the
  # pooled SD itself must be representable for the objective to be finite.
  values <- c(fn, scale)
  powers <- pmin(floor(log2(values)), 1023)
  multiplicity <- c(1, rep(power, length(scale)))
  fraction <- prod((values / 2^powers)^multiplicity)
  exponent <- sum(powers * multiplicity)
  base <- max(-1022, min(1023, exponent))
  (fraction * 2^(exponent - base)) * 2^base
}

#' Dispatch to the C++ solver for a given method
#' @keywords internal
#' @noRd
reg_solver <- function(method, ..., pooled_ratio = 0,
                       save_weights = c(1, 1, 1)) {
  switch(method,
    sir   = sir_solver(..., pooled_ratio = pooled_ratio),
    save  = save_solver(..., identity_weight = save_weights[1L],
                        covariance_weight = save_weights[2L], moment_scale = save_weights[3L]),
    phd   = phd_solver(...),
    seff  = seff_solver(...),
    local = local_solver(...),
    stop("Unknown method: ", method)
  )
}

#' Evaluate the estimating-equation objective at a candidate B (for init
#' selection).
#' @keywords internal
#' @noRd
reg_init <- function(method, ..., pooled_ratio = 0,
                     save_weights = c(1, 1, 1)) {
  switch(method,
    sir   = sir_init(..., pooled_ratio = pooled_ratio),
    save  = save_init(..., identity_weight = save_weights[1L],
                      covariance_weight = save_weights[2L], moment_scale = save_weights[3L]),
    phd   = phd_init(...),
    seff  = seff_init(...),
    local = local_f(...),
    stop("Unknown method: ", method)
  )
}

#' Pick the best starting value for B
#'
#' Evaluates the estimating-equation objective at SIR, SAVE, and PHD initial
#' B matrices and returns the one with the smallest value.  Used by all
#' methods except \code{method = "local"} with \code{ndr = 1}, which has
#' its own multi-start strategy with NW-\eqn{R^2} selection.
#'
#' For two-index models (\code{ndr = 2}), a sequential candidate
#' combining the leading direction from SIR and SAVE is added to help span
#' the joint subspace.
#'
#' @keywords internal
#' @noRd
initB <- function(x, y, ndr, bw, method, ncore) {
  candidates <- reg_candidates(x, y, ndr)
  select_reg_candidate(candidates, x, y, bw, method, ncore, pooled_ratio = 0)
}

reg_candidates <- function(x, y, ndr) {
  B1 <- gramSchmidt(dr(y ~ x, method = "sir")$evectors[, 1:ndr, drop = FALSE])$Q
  B2 <- gramSchmidt(dr(y ~ x, method = "save")$evectors[, 1:ndr, drop = FALSE])$Q
  B3 <- gramSchmidt(dr(y ~ x, method = "phd")$evectors[, 1:ndr, drop = FALSE])$Q

  candidates <- list(B1, B2, B3)

  # For two-index models, add a sequential candidate
  # combining the leading direction from SIR and SAVE.  Each sliced
  # method captures different structure (mean vs. variance), so the
  # joint subspace is often better spanned by mixed directions.
  if (ndr == 2L) {
    B_seq <- cbind(
      dr(y ~ x, method = "sir")$evectors[, 1],
      dr(y ~ x, method = "save")$evectors[, 1]
    )
    if (all(is.finite(B_seq)) && qr(B_seq)$rank == ndr)
      candidates[[4]] <- gramSchmidt(B_seq)$Q
  }
  candidates
}

select_reg_candidate <- function(candidates, x, y, bw, method, ncore, pooled_ratio,
                                 save_weights = c(1, 1, 1)) {
  values <- vapply(candidates, function(Bc) {
    reg_init(method, Bc, x, y, bw, ncore, pooled_ratio = pooled_ratio,
             save_weights = save_weights)
  }, numeric(1))

  candidates[[which.min(values)]]
}
