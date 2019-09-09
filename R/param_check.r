#' @title param_check
#' @name param_check
#' @description Check parameters
#' @param control control parameters
#' @keywords internal


control.check <- function(control) {
  if (!is.list(control)) {
    stop("control must be a list of tuning parameters")
  }

  if (is.null(control$rho)) {
    control$rho <- 1e-4
  } else if (control$rho < 0) control$rho <- 1e-4

  if (is.null(control$eta)) {
    control$eta <- 0.2
  } else if (control$eta < 0) control$eta <- 0.2

  if (is.null(control$gamma)) {
    control$gamma <- 0.85
  } else if (control$gamma < 0) control$gamma <- 0.85

  if (is.null(control$tau)) {
    control$tau <- 1e-3
  } else if (control$tau < 0) control$tau <- 1e-3

  if (is.null(control$epsilon)) {
    control$epsilon <- 1e-6
  } else if (control$epsilon < 0) control$epsilon <- 1e-6

  if (is.null(control$btol)) {
    control$btol <- 1e-6
  } else if (control$btol < 0) control$btol <- 1e-6

  if (is.null(control$ftol)) {
    control$ftol <- 1e-6
  } else if (control$ftol < 0) control$ftol <- 1e-6

  if (is.null(control$gtol)) {
    control$gtol <- 1e-6
  } else if (control$gtol < 0) control$gtol <- 1e-6

  return(control)
}
