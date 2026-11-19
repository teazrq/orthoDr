#ifndef ORTHODR_SOLVER_NUMERICS_H
#define ORTHODR_SOLVER_NUMERICS_H

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

namespace orthoDr {

class CayleyCurve {
  bool dense;
  arma::mat H, HB, U, VU, VB;

public:
  CayleyCurve(const arma::mat& B, const arma::mat& G, bool use_dense)
    : dense(use_dense) {
    if (dense) {
      const arma::mat GB = G * B.t();
      H = 0.5 * (GB - GB.t());
      HB = H * B;
    } else {
      U = arma::join_rows(G, B);
      const arma::mat V = arma::join_rows(B, -G);
      VU = V.t() * U;
      VB = V.t() * B;
    }
  }

  bool trial(const arma::mat& B, double step, arma::mat& result) const {
    if (dense) {
      const arma::mat lhs = arma::eye(H.n_rows, H.n_cols) + step * H;
      const arma::mat rhs = B - step * HB;
      return lhs.is_finite() && rhs.is_finite() &&
        arma::solve(result, lhs, rhs, arma::solve_opts::no_approx);
    }
    const arma::mat lhs = arma::eye(VU.n_rows, VU.n_cols) + 0.5 * step * VU;
    arma::mat coefficients;
    if (!lhs.is_finite() || !arma::solve(coefficients, lhs, VB, arma::solve_opts::no_approx))
      return false;
    result = B - U * (step * coefficients);
    return true;
  }
};

inline double feasibility_error(const arma::mat& B) {
  return arma::norm(B.t() * B - arma::eye(B.n_cols, B.n_cols), "fro");
}

constexpr double feasibility_tolerance = 1e-10;

template<class Objective>
bool backtrack(arma::mat& B, double& F, double& step, double eta,
               double reference, double slope, const CayleyCurve& curve,
               const Objective& objective) {
  for (int trial = 0; trial < 30; ++trial) {
    if (!std::isfinite(step) || step <= 0) break;
    const double threshold = reference - step * slope;
    arma::mat candidate;
    if (std::isfinite(threshold) && curve.trial(B, step, candidate) &&
        candidate.is_finite() && feasibility_error(candidate) <= feasibility_tolerance &&
        !arma::approx_equal(candidate, B, "absdiff", 0.0)) {
      const double value = objective(candidate);
      if (std::isfinite(value) && value <= threshold) {
        B = std::move(candidate);
        F = value;
        return true;
      }
    }
    const double smaller = eta * step;
    if (smaller >= step) break;
    step = smaller;
  }
  return false;
}

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
