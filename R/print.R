#' @title Print an orthoDr Fit or Prediction Summary
#' @description Prints the estimated direction matrix for a fitted model,
#'   or identifies the returned quantities for a prediction object.
#' @export
#' @keywords internal
#' @param x An orthoDr fitted model or prediction object.
#' @param ... Additional arguments, currently unused.
#' @return Invisible NULL. Called for its side-effect (printing to console).
#' @examples
#' # generate some survival data
#' N = 100; P = 4; dataX = matrix(rnorm(N*P), N, P)
#' Y = exp(-1 + dataX[,1] + rnorm(N))
#' Censor = rbinom(N, 1, 0.8)
#'
#' # fit the model
#' orthoDr_surv(dataX, Y, Censor, ndr = 1, method = "dm")

print.orthoDr <- function(x, ...)
{
  if (inherits(x, "fit"))
  {
    cat(paste("Subspace for", class(x)[3], "model using", x$method, "approach:\n"))
    print(x$B)
  }

  if (inherits(x, "predict"))
  {

    if (inherits(x, "reg"))
      cat(paste("Prediction for orthoDr regression: mean prediction"))

    if (inherits(x, "surv"))
    {
      cat(paste("Prediction for orthoDr Survival:", ncol(x$surv), "testing subjects at", length(x$timepoints), "time points\n"))
      cat("See 'surv' and 'timepoints'.")
    }

    if (inherits(x, "pdose"))
      cat(paste("Prediction for orthoDr personalized dose: predicted doses"))
  }
}







