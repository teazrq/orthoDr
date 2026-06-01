#' Semiparametric Dimension Reduction for Regression
#'
#' @description
#' Fits a semiparametric sufficient dimension reduction model for continuous
#' outcomes using the methods of Ma & Zhu (2012, 2013). The direction matrix
#' \code{B} is estimated via orthogonality-constrained optimization on the
#' Stiefel manifold, ensuring \eqn{B^\top B = I}.
#'
#' Kernel density estimation operates on the scaled projections \eqn{\beta^\top x},
#' so the procedure is robust to the scale of \code{x}. However, raw feature
#' values enter the Nadaraya–Watson estimator unscaled — users with features of
#' very different variance should standardise \code{x} beforehand. The outcome
#' \code{y} is internally scaled.
#'
#' @name orthoDr_reg
#' @export
#'
#' @param x A numeric matrix of predictors (n observations \eqn{\times} p
#'   features). Not automatically scaled.
#' @param y A numeric vector of continuous outcomes (length n).
#' @param method The dimension reduction method. One of \code{"sir"},
#'   \code{"save"}, \code{"phd"}, \code{"local"}, or \code{"seff"}.
#'   See Details for method-specific behaviour.
#' @param ndr Integer; number of directions to estimate (default 2). Must be
#'   \eqn{\ge 1} and \eqn{\le} \code{ncol(x)}. Values > 4 produce a warning.
#' @param B.initial Initial \eqn{B} matrix (\eqn{p \times ndr}), or \code{NULL}
#'   (default). When \code{NULL}, \code{initB()} selects the best among SIR,
#'   SAVE, and PHD starters by minimising the squared norm of the estimating
#'   equation.  **Exception:** for \code{method = "local"} with \code{ndr = 1},
#'   a multi-start strategy runs LOCAL from all three inits and selects the
#'   result that maximises the nonparametric \eqn{R^2} of the Nadaraya–Watson
#'   regression of \eqn{Y} on \eqn{B^T X}.
#' @param bw X-kernel bandwidth. \code{NULL} (default) uses Silverman's rule
#'   of thumb. Internal Y-kernel bandwidths are derived from the sample size
#'   and structural dimension — they are not user-tunable.
#' @param keep.data Logical; if \code{TRUE}, the original data is stored in the
#'   returned object for use by \code{\link{predict}}. Default \code{FALSE}.
#' @param control A named list of optimizer tolerances. Recognised entries and
#'   their defaults:
#'   \describe{
#'     \item{\code{rho}}{1e-4 — augmented Lagrangian penalty}
#'     \item{\code{eta}}{0.2 — step-size decay factor}
#'     \item{\code{gamma}}{0.85 — linesearch contraction}
#'     \item{\code{tau}}{1e-3 — initial step-size adjustment}
#'     \item{\code{epsilon}}{1e-6 — finite-difference perturbation}
#'     \item{\code{btol}}{1e-6 — parameter-change tolerance}
#'     \item{\code{ftol}}{1e-6 — objective-change tolerance}
#'     \item{\code{gtol}}{1e-6 — gradient-norm tolerance}
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
#'     \item{PHD}{Semiparametric principal Hessian directions.}
#'     \item{Local}{Locally efficient estimator (Ma & Zhu 2013,
#'       Section~3.1) with a normal posited model for the conditional
#'       density.  This yields a mean-regression estimator — it captures
#'       dependence of \eqn{Y} on \eqn{\beta^\top x} through the conditional
#'       mean only, not higher moments.  For \eqn{d = 1} with default
#'       initialisation, a multi-start strategy runs LOCAL from SIR, SAVE,
#'       and PHD starters and selects the result that maximises the
#'       nonparametric \eqn{R^2} of the Nadaraya–Watson regression of
#'       \eqn{Y} on \eqn{\beta^\top x}.}
#'     \item{SEFF}{Semiparametric efficient estimator (Ma & Zhu 2013,
#'       Section 3.2). Uses a two-step procedure: (1) obtain a root-\eqn{n}
#'       consistent initial \eqn{\tilde\beta} via SAVE; (2) estimate the
#'       nonparametric nuisance components (conditional density
#'       \eqn{\hat\eta_2} and its derivative, plus \eqn{\hat{E}[x \mid
#'       \beta^\top x]}) once at \eqn{\tilde\beta}; (3) optimise the
#'       efficient score with those estimates held fixed.}
#'   }
#' }
#'
#' @return An object of class \code{c("orthoDr", "fit", "reg")}, a list with
#'   components:
#'   \describe{
#'     \item{\code{B}}{Estimated \eqn{p \times ndr} direction matrix with
#'       orthonormal columns.}
#'     \item{\code{fn}}{Final objective value.}
#'     \item{\code{itr}}{Number of iterations performed.}
#'     \item{\code{converge}}{Convergence code (0 = success).}
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
#' # Mean model — SIR
#' Y <- -1 + X[, 1] + rnorm(N)
#' fit_sir <- orthoDr_reg(X, Y, ndr = 1, method = "sir")
#' fit_sir$B
#'
#' # Variance model — PHD
#' Y <- -1 + X[, 1]^2 + rnorm(N)
#' fit_phd <- orthoDr_reg(X, Y, ndr = 1, method = "phd")
#' fit_phd$B
#'
#' # Efficient estimator (SEFF) — works on both mean and variance models
#' fit_seff  <- orthoDr_reg(X, Y, ndr = 1, method = "seff")
#' fit_seff$B
#'
#' # Custom initial + prediction
#' B0 <- matrix(c(1, 0, 0, 0), 4, 1)
#' fit <- orthoDr_reg(X, Y, ndr = 1, B.initial = B0, keep.data = TRUE)
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

  if (is.null(bw))
    bw <- silverman(ndr, N)

  # --- data preparation ---
  # y is internally scaled per Ma & Zhu convention
  Y <- as.matrix(scale(y) / N^(-1 / (ndr + 5)) / sqrt(2))

  # x is NOT scaled — kernel operates on scaled βᵀx, not raw x.
  # Users with features of disparate scale should standardise beforehand.
  X <- x

  # --- initial value ---
  use_default_init <- is.null(B.initial)
  if (use_default_init) {
    B.initial <- initB(X, Y, ndr, bw, method, ncore)
  } else {
    if (!is.matrix(B.initial))
      stop("B.initial must be a matrix")
    if (ncol(x) != nrow(B.initial) || ndr != ncol(B.initial))
      stop("Dimension of B.initial is not correct")
    B.initial <- gramSchmidt(B.initial)$Q
  }

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

    B_sir  <- gramSchmidt(dr(y ~ x, method = "sir")$evectors[, 1:ndr, drop = FALSE])$Q
    B_save <- gramSchmidt(dr(y ~ x, method = "save")$evectors[, 1:ndr, drop = FALSE])$Q
    B_phd  <- gramSchmidt(dr(y ~ x, method = "phd")$evectors[, 1:ndr, drop = FALSE])$Q

    fit_sir <- reg_solver(
      "local", B_sir, X, Y, bw,
      control$rho, control$eta, control$gamma,
      control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol,
      maxitr, 0L, ncore
    )
    fit_save <- reg_solver(
      "local", B_save, X, Y, bw,
      control$rho, control$eta, control$gamma,
      control$tau, control$epsilon,
      control$btol, control$ftol, control$gtol,
      maxitr, 0L, ncore
    )
    fit_phd <- reg_solver(
      "local", B_phd, X, Y, bw,
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

  fit <- reg_solver(
    method, B.initial, X, Y, bw,
    control$rho, control$eta, control$gamma,
    control$tau, control$epsilon,
    control$btol, control$ftol, control$gtol,
    maxitr, verbose, ncore
  )

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

#' Dispatch to the C++ solver for a given method
#' @keywords internal
#' @noRd
reg_solver <- function(method, ...) {
  switch(method,
    sir   = sir_solver(...),
    save  = save_solver(...),
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
reg_init <- function(method, ...) {
  switch(method,
    sir   = sir_init(...),
    save  = save_init(...),
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
#' For multi-index models (\code{ndr > 1}), a sequential candidate
#' combining the leading direction from SIR and SAVE is added to help span
#' the joint subspace.
#'
#' @keywords internal
#' @noRd
initB <- function(x, y, ndr, bw, method, ncore) {

  B1 <- gramSchmidt(dr(y ~ x, method = "sir")$evectors[, 1:ndr, drop = FALSE])$Q
  B2 <- gramSchmidt(dr(y ~ x, method = "save")$evectors[, 1:ndr, drop = FALSE])$Q
  B3 <- gramSchmidt(dr(y ~ x, method = "phd")$evectors[, 1:ndr, drop = FALSE])$Q

  candidates <- list(B1, B2, B3)

  # For multi-index models (ndr > 1), add a sequential candidate
  # combining the leading direction from SIR and SAVE.  Each sliced
  # method captures different structure (mean vs. variance), so the
  # joint subspace is often better spanned by mixed directions.
  if (ndr > 1) {
    B_seq <- gramSchmidt(cbind(
      dr(y ~ x, method = "sir")$evectors[, 1],
      dr(y ~ x, method = "save")$evectors[, 1]
    ))$Q
    candidates[[4]] <- B_seq
  }

  values <- vapply(candidates, function(Bc) {
    reg_init(method, Bc, x, y, bw, ncore)
  }, numeric(1))

  candidates[[which.min(values)]]
}
