#' @title orthoDr: Semi-Parametric Dimension Reduction Models Using Orthogonality Constrained Optimization
#' @description Provides semi-parametric dimension reduction methods using
#'   orthogonality constrained optimization. Includes regression, survival,
#'   and personalized dose models.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

# The following blocks are used by roxygen2 to generate NAMESPACE imports.
# They are placed here so all package-level imports are centralized.

#' @importFrom Rcpp evalCpp
#' @importFrom pracma gramSchmidt repmat
#' @importFrom stats pnorm cancor cov dist weighted.mean sd quantile predict
#' @importFrom survival survfit Surv
#' @importFrom plot3D mesh
#' @importFrom rgl surface3d axis3d mtext3d box3d
#' @importFrom grDevices rainbow
#' @importFrom graphics legend par plot plot.new
#' @importFrom dr dr
#' @importFrom MASS ginv mvrnorm
NULL

#' @useDynLib orthoDr, .registration = TRUE
NULL
