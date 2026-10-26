#ifndef ORTHODR_SOLVER_NUMERICS_H
#define ORTHODR_SOLVER_NUMERICS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

namespace orthoDr {

inline void check_state(double F, const arma::mat& G) {
  if (!std::isfinite(F)) Rcpp::stop("Objective must be finite at every evaluation");
  if (!G.is_finite()) Rcpp::stop("Gradient must be finite at every evaluation");
}

// Preserve the solver's bounds without rejecting small, well-defined ratios.
inline double bb_step(double SS, double SY, double YY, bool even,
                      double previous, double upper, bool& valid) {
  const double numerator = even ? SS : SY;
  const double denominator = even ? SY : YY;
  valid = std::isfinite(numerator) && numerator > 0 &&
    std::isfinite(denominator) && denominator > 0;
  double candidate = valid ? numerator / denominator : previous;
  valid = valid && std::isfinite(candidate) && candidate > 0;
  if (!valid) candidate = previous;
  if (!std::isfinite(candidate) || candidate <= 0)
    Rcpp::stop("Optimizer step must be positive and finite");
  return std::max(1e-20, std::min(candidate, upper));
}

} // namespace orthoDr

#endif
