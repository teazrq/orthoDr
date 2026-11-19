#' Control parameter checking
#' @param control control parameters
#' @name control.check
#' @keywords internal
#' @noRd
control.check <- function(control)
{
  if (!is.list(control))
    stop("control must be a list of tuning parameters")

  defaults <- list(
    rho     = 1e-4,
    eta     = 0.2,
    gamma   = 0.85,
    tau     = 1e-3,
    epsilon = 1e-6,
    btol    = 1e-6,
    ftol    = 1e-6,
    gtol    = 1e-6
  )

  for (name in names(defaults)) {
    val <- control[[name]]
    if (is.null(val) || !is.finite(val) || val <= 0)
      control[[name]] <- defaults[[name]]
  }

  return(control)
}

orthonormalize_start <- function(B)
{
  Q <- gramSchmidt(B)$Q
  # Keep the original rank check and refine starts that miss the solver tolerance.
  if (norm(crossprod(Q) - diag(ncol(Q)), "F") > 1e-10)
    Q <- gramSchmidt(Q)$Q
  Q
}

check_survival_censor <- function(censor, n, retained = FALSE)
{
  advice <- if (retained) " Recode the event indicators and refit the model." else ""
  coding <- "Use 0 (censored) or 1 (event)."
  if (!(is.numeric(censor) || is.logical(censor)) ||
      !(is.null(dim(censor)) || (is.matrix(censor) && ncol(censor) == 1L)))
    stop("censor must be a numeric or logical vector or a one-column matrix. ",
         coding, advice, call. = FALSE)
  if (length(censor) != n)
    stop("Number of observations do not match. censor must have one value per observation. ",
         coding, advice, call. = FALSE)

  invalid <- which(!censor %in% c(0, 1))
  if (length(invalid)) {
    shown <- invalid[seq_len(min(length(invalid), 8L))]
    examples <- paste0(shown, " (", sprintf("%.17g", censor[shown]), ")", collapse = ", ")
    remaining <- length(invalid) - length(shown)
    more <- if (remaining) paste0(" and ", remaining, " more") else ""
    stop("Invalid censor at 1-based observation rows ", examples, more, ". ",
         coding, advice, call. = FALSE)
  }
}

check_survival_time <- function(y, n, name = "y", retained = FALSE)
{
  advice <- if (retained) " Correct the original follow-up times and refit the survival model." else ""
  if (!((is.numeric(y) && !is.complex(y)) || is.logical(y) ||
        inherits(y, c("Date", "POSIXt", "difftime"))) ||
      !(is.null(dim(y)) || (is.matrix(y) && ncol(y) == 1L)))
    stop(name, " must be a real numeric or time vector or a one-column matrix.",
         advice, call. = FALSE)
  if (length(y) != n)
    stop("Number of observations do not match. ", name,
         " must have one value per observation.", advice, call. = FALSE)

  invalid <- which(!is.finite(y))
  if (length(invalid)) {
    shown <- invalid[seq_len(min(length(invalid), 8L))]
    remaining <- length(invalid) - length(shown)
    more <- if (remaining) paste0(" and ", remaining, " more") else ""
    stop(name, " must contain finite times. Invalid at 1-based observation rows ",
         paste(shown, collapse = ", "), more, ".", advice, call. = FALSE)
  }
}

survival_coordinates <- function(x)
{
  if (!is.numeric(x) || is.complex(x) || any(!is.finite(x)) ||
      nrow(x) < 2L || ncol(x) < 1L)
    stop("Survival fitting requires a finite numeric predictor matrix with at least two rows.")
  delta <- sweep(x, 2L, as.double(x[1L, ]), "-")
  span <- apply(abs(delta), 2L, max)
  if (any(span == 0))
    stop("Constant predictor at column(s) ", paste(which(span == 0), collapse = ", "),
         ". Each predictor must have positive sample variance.")
  if (any(!is.finite(span))) stop("Predictor differences exceed the numeric range. Rescale x.")
  relative <- sweep(delta, 2L, span, "/")
  centered <- sweep(relative, 2L, colMeans(relative), "-")
  sd.relative <- sqrt(colSums(centered^2) / (nrow(x) - 1L))
  scale <- span * sd.relative
  if (any(!is.finite(scale) | scale == 0)) stop("Predictor scales exceed the numeric range. Rescale x.")
  list(x = sweep(centered, 2L, sd.relative, "/"), scale = scale)
}

survival_design_svd <- function(x, nv = 0L)
{
  s <- svd(x, nu = 0L, nv = nv)
  tol <- max(dim(x)) * .Machine$double.eps * s$d[1L]
  if (length(s$d) < ncol(x) || any(s$d <= tol))
    stop("Survival fitting requires full centered predictor rank (a positive definite feature covariance). Remove redundant predictors.")
  s
}

raw_predictor_directions <- function(b, scale)
{
  message <- "Cannot express the fitted direction in the supplied predictor units. Rescale x."
  for (j in seq_len(ncol(b))) {
    v <- b[, j]
    if (any(!is.finite(v)) || max(abs(v)) == 0) stop(message)
    raw <- v / scale
    if (any(!is.finite(raw)) || max(abs(raw)) == 0) {
      mag <- log(abs(v)) - log(scale)
      raw <- sign(v) * exp(mag - max(mag))
    }
    raw <- raw / max(abs(raw))
    raw <- raw / sqrt(sum(raw^2))
    # Check the represented unit vector back in the working coordinates.
    mag <- log(abs(raw)) + log(scale)
    restored <- sign(raw) * exp(mag - max(mag))
    if (max(abs(restored - v / max(abs(v)))) > sqrt(.Machine$double.eps)) stop(message)
    b[, j] <- raw
  }
  b
}

check_survival_projection <- function(x, b)
{
  s <- svd(x %*% b, nu = 0L, nv = 0L)$d
  tol <- max(dim(x)) * .Machine$double.eps * norm(x, "F")
  if (length(s) < ncol(b) || any(s <= tol))
    stop("The survival directions must give linearly independent, varying projected predictors.")
}
