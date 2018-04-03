#' @title print.orthoDr
#' @name print
#' @description Print a orthoDr object
#' @keywords internal
#' @param x A fitted orthoDr object
#' @param ... ...
#' @examples
#' # generate some survival data
#' N = 100; P = 4; dataX = matrix(rnorm(N*P), N, P)
#' Y = exp(-1 + dataX[,1] + rnorm(N))
#' Censor = rbinom(N, 1, 0.8)

#' # fit the model
#' orthoDr_surv(dataX, Y, Censor, ndr = 1, method = "dm")

print.orthoDr<- function(x, ...)
{
  if (class(x)[2] == "fit")
  {
    cat(paste("Subspace for", class(x)[3], "model using", x$method, "approach:\n"))
    print(x$B)
  }

  if (class(x)[2] == "predict")
  {
    if (class(x)[3] == "survival")
    {
      cat(paste("Prediction for orthoDr Survival:", ncol(x$surv), "testing subjects at", length(x$timepoints), "time points\n"))
      cat("See 'surv' and 'timepoints'.")
    }
    if (class(x)[3] == "pt")
      cat(paste("Prediction for orthoDr personalized treatment: best treatment dose and reward prediction"))

    if (class(x)[3] == "regression")
      cat(paste("Prediction for orthoDr regression: mean prediction"))
  }
}







