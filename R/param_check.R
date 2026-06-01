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
