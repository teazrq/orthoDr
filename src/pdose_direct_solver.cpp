//    ----------------------------------------------------------------
//
//    Orthogonality Constrained Optimization for Dimension Reduction
//    (orthoDr)
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 2
//    of the License, or (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public
//    License along with this program; if not, write to the Free
//    Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
//    Boston, MA  02110-1301, USA.
//
//    ----------------------------------------------------------------

#include <RcppArmadillo.h>
#include "utilities.h"
#include "orthoDr_pdose.h"
#include "solver_numerics.h"
#include <limits>


using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

struct DirectKernel {
  arma::mat projected;
  arma::mat centered;
  arma::rowvec scales;
  arma::mat log_kernel;
  arma::mat kernel;
  arma::rowvec mean;
};

struct DirectPolicy {
  bool intercept;
  double lower, upper;

  DirectPolicy() : intercept(false), lower(0), upper(0) {}
  double project(const double alpha) const {
    return std::max(lower, std::min(upper, alpha));
  }
};

DirectPolicy direct_policy(const bool intercept,
                            const Rcpp::Nullable<Rcpp::NumericVector>& bounds,
                            const arma::vec& A) {
  DirectPolicy out;
  out.intercept = intercept;
  if (!intercept) {
    if (bounds.isNotNull()) Rcpp::stop("intercept.bounds requires intercept = TRUE");
    return out;
  }
  if (bounds.isNotNull()) {
    const Rcpp::NumericVector limits(bounds.get());
    if (limits.size() != 2 || !std::isfinite(limits[0]) ||
        !std::isfinite(limits[1]) || limits[0] >= limits[1])
      Rcpp::stop("intercept.bounds must contain two finite, strictly ordered values");
    out.lower = limits[0];
    out.upper = limits[1];
  } else {
    if (A.n_elem == 0 || !A.is_finite()) Rcpp::stop("Observed doses must be finite");
    out.lower = A.min();
    out.upper = A.max();
    if (out.lower >= out.upper)
      Rcpp::stop("Default intercept.bounds requires varying observed doses");
  }
  return out;
}

struct DirectEvaluation {
  bool valid;
  std::string message;
  double objective;
  double value;
  double penalty;
  double norm_B;
  double norm_W;
  bool intercept;
  double alpha, gradient_alpha, residual_alpha, norm_policy;
  arma::vec dose;
  arma::mat features;
  arma::mat gradient_B;
  arma::vec gradient_W;
  arma::vec curvature_dose;
  DirectKernel geometry;

  DirectEvaluation() : valid(false), objective(NA_REAL), value(NA_REAL),
    penalty(NA_REAL), norm_B(NA_REAL), norm_W(NA_REAL), intercept(false),
    alpha(0), gradient_alpha(NA_REAL), residual_alpha(NA_REAL), norm_policy(NA_REAL) {}
};

bool direct_kernel(const arma::mat& B, const arma::mat& X,
                   const double bw, DirectKernel& out, std::string& message) {
  if (X.n_rows < 2 || X.n_cols == 0 || B.n_rows != X.n_cols ||
      B.n_cols == 0 || B.n_cols > B.n_rows || !X.is_finite() ||
      !B.is_finite() || !std::isfinite(bw) || bw <= 0) {
    message = "Direct-dose data, basis dimensions, and bandwidth must be valid and finite";
    return false;
  }
  out.projected = X * B;
  if (!out.projected.is_finite()) {
    message = "Projected covariates must be finite";
    return false;
  }
  out.centered = out.projected.each_row() - arma::mean(out.projected, 0);
  out.scales = arma::stddev(out.projected, 0, 0);
  if (!out.centered.is_finite() || !out.scales.is_finite() ||
      arma::any(out.scales <= 0)) {
    message = "Projected covariates must have positive finite standard deviations";
    return false;
  }
  const arma::uword n = X.n_rows;
  out.log_kernel.zeros(n, n);
  for (arma::uword k = 0; k < B.n_cols; ++k) {
    Rcpp::checkUserInterrupt();
    const double scale = bw * out.scales(k);
    if (!std::isfinite(scale) || scale <= 0) {
      message = "Projected bandwidth scales must be positive and finite";
      return false;
    }
    arma::vec z = out.centered.col(k) / scale;
    arma::mat delta = arma::repmat(z, 1, n) - arma::repmat(z.t(), n, 1);
    out.log_kernel -= 0.5 * arma::square(delta);
    if (!out.log_kernel.is_finite()) {
      message = "Covariate log-kernel arithmetic is nonfinite";
      return false;
    }
  }
  out.kernel = arma::exp(out.log_kernel);
  out.mean = arma::mean(out.kernel, 0);
  return out.kernel.is_finite();
}

DirectEvaluation direct_evaluate(const arma::mat& B, const arma::mat& X,
                                 const arma::vec& A, const arma::vec& R,
                                 const double bw, const arma::vec& W,
                                 const double lambda_reward,
                                 const bool gradients = true,
                                 const DirectPolicy& policy = DirectPolicy(),
                                 const double alpha = 0) {
  DirectEvaluation out;
  out.intercept = policy.intercept;
  out.alpha = alpha;
  const arma::uword n = X.n_rows;
  if (A.n_elem != n || R.n_elem != n || W.n_elem != n ||
      !A.is_finite() || !R.is_finite() || !W.is_finite() ||
      !std::isfinite(lambda_reward) || lambda_reward < 0) {
    out.message = "Direct-dose outcomes and coefficients must be finite with matching lengths, and lambda.reward must be nonnegative";
    return out;
  }
  if (policy.intercept && (lambda_reward <= 0 || !std::isfinite(alpha) ||
      alpha < policy.lower || alpha > policy.upper)) {
    out.message = "Intercept mode requires lambda.reward > 0 and a finite intercept within its bounds";
    return out;
  }
  if (!direct_kernel(B, X, bw, out.geometry, out.message)) return out;
  const double dose_scale = arma::stddev(A) * bw;
  if (!std::isfinite(dose_scale) || dose_scale <= 0) {
    out.message = "Observed doses must have a positive finite bandwidth scale";
    return out;
  }
  out.features = out.geometry.kernel;
  if (policy.intercept) out.features.each_row() -= out.geometry.mean;
  out.dose = out.features * W;
  if (policy.intercept) out.dose += alpha;
  if (!out.dose.is_finite()) {
    out.message = "Policy doses are nonfinite";
    return out;
  }
  arma::mat log_sensitivity;
  arma::vec gradient_dose;
  if (gradients) {
    log_sensitivity.set_size(n, n);
    gradient_dose.set_size(n);
    out.curvature_dose.set_size(n);
  }
  out.value = 0;
  for (arma::uword i = 0; i < n; ++i) {
    if (i % 32 == 0) Rcpp::checkUserInterrupt();
    const arma::vec dose_delta = (A - out.dose(i)) / dose_scale;
    const arma::vec log_weight = out.geometry.log_kernel.col(i) - arma::square(dose_delta);
    if (!log_weight.is_finite()) {
      out.message = "Reward log-weight arithmetic is nonfinite";
      return out;
    }
    arma::vec probability = arma::exp(log_weight - log_weight.max());
    const double normalizer = arma::accu(probability);
    if (!std::isfinite(normalizer) || normalizer <= 0) {
      out.message = "Reward weights cannot be normalized";
      return out;
    }
    probability /= normalizer;
    // Centering gives exactly zero reward derivatives for constant rewards.
    const double reward = R(0) + arma::dot(probability, R - R(0));
    out.value += reward / static_cast<double>(n);
    if (gradients) {
      const arma::vec alpha = probability % (R - reward);
      log_sensitivity.col(i) = -alpha / static_cast<double>(n);
      gradient_dose(i) = -2.0 * arma::dot(alpha, dose_delta) /
        (static_cast<double>(n) * dose_scale);
      const arma::vec centered_dose = dose_delta - arma::dot(probability, dose_delta);
      // Optional curvature may overflow even when the objective and first
      // derivatives are finite. Such a state still permits a gradient step.
      out.curvature_dose(i) =
        (-4.0 * arma::dot(alpha, arma::square(centered_dose)) /
         (static_cast<double>(n) * dose_scale)) / dose_scale;
    }
  }
  out.penalty = lambda_reward == 0 ? 0 : lambda_reward * arma::dot(W, W);
  out.objective = -out.value + out.penalty;
  if (!std::isfinite(out.value) || !std::isfinite(out.penalty) ||
      !std::isfinite(out.objective)) {
    out.message = "Direct-dose objective components are nonfinite";
    return out;
  }
  if (gradients) {
    out.gradient_W = out.features.t() * gradient_dose +
      2.0 * (lambda_reward * W);
    out.gradient_alpha = policy.intercept ? arma::accu(gradient_dose) : 0;
    out.residual_alpha = policy.intercept ?
      alpha - policy.project(alpha - out.gradient_alpha) : 0;
    arma::vec policy_gradient = gradient_dose;
    if (policy.intercept) policy_gradient -= arma::mean(gradient_dose);
    // Each column i has dose d_i = sum_j K(j,i) W_j. Both kernel
    // endpoints and the sample standard deviation move with B.
    const arma::mat sensitivity = log_sensitivity +
      out.geometry.kernel % (W * policy_gradient.t());
    const arma::mat symmetric_sensitivity = sensitivity + sensitivity.t();
    arma::mat gradient_projected(n, B.n_cols);
    for (arma::uword k = 0; k < B.n_cols; ++k) {
      Rcpp::checkUserInterrupt();
      const double scale = out.geometry.scales(k);
      const arma::vec z = out.geometry.centered.col(k) / (bw * scale);
      const arma::mat delta = arma::repmat(z, 1, n) - arma::repmat(z.t(), n, 1);
      const double scale_sensitivity = arma::accu(sensitivity % arma::square(delta));
      gradient_projected.col(k) =
        -arma::sum(symmetric_sensitivity % delta, 1) / (bw * scale) +
        scale_sensitivity * (out.geometry.centered.col(k) / scale) /
        ((static_cast<double>(n) - 1.0) * scale);
    }
    out.gradient_B = X.t() * gradient_projected;
    const arma::mat BtG = B.t() * out.gradient_B;
    const arma::mat projected_gradient = out.gradient_B - B * (0.5 * (BtG + BtG.t()));
    out.norm_B = arma::norm(projected_gradient, "fro");
    out.norm_W = arma::norm(out.gradient_W, 2);
    out.norm_policy = std::hypot(out.norm_W, out.residual_alpha);
    if (!out.gradient_B.is_finite() || !out.gradient_W.is_finite() ||
        !std::isfinite(out.norm_B) || !std::isfinite(out.norm_policy) ||
        !std::isfinite(out.gradient_alpha)) {
      out.message = "Direct-dose derivatives are nonfinite";
      return out;
    }
  }
  out.valid = true;
  return out;
}

bool direct_hessian_W(const DirectEvaluation& state, const double lambda_reward,
                       arma::mat& hessian) {
  if (!state.valid || state.curvature_dose.n_elem != state.geometry.kernel.n_rows ||
      !state.curvature_dose.is_finite()) {
    hessian.set_size(0, 0);
    return false;
  }
  arma::mat weighted_kernel = state.features;
  weighted_kernel.each_col() %= state.curvature_dose;
  hessian = state.features.t() * weighted_kernel;
  hessian.diag() += 2.0 * lambda_reward;
  // The exact expression is symmetric. Copy one triangle rather than
  // introducing rounding asymmetry into the Cholesky factorization.
  hessian = arma::symmatu(hessian);
  return hessian.is_finite();
}

arma::vec direct_direction_W(const DirectEvaluation& state,
                              const double lambda_reward,
                              std::string& type) {
  arma::vec direction = -state.gradient_W;
  type = "gradient";
  arma::mat hessian, factor;
  arma::vec intermediate, solution;
  if (direct_hessian_W(state, lambda_reward, hessian) &&
      arma::chol(factor, hessian) &&
      arma::solve(intermediate, arma::trimatl(factor.t()), state.gradient_W,
                   arma::solve_opts::no_approx) &&
      arma::solve(solution, arma::trimatu(factor), intermediate,
                   arma::solve_opts::no_approx)) {
    const double slope = -arma::dot(state.gradient_W, solution);
    if (solution.is_finite() && std::isfinite(slope) && slope < 0) {
      direction = -solution;
      type = "newton";
    }
  }
  return direction;
}

bool direct_hessian_policy(const DirectEvaluation& state,
                            const double lambda_reward, arma::mat& hessian) {
  arma::mat coefficient;
  if (!direct_hessian_W(state, lambda_reward, coefficient)) {
    hessian.set_size(0, 0);
    return false;
  }
  const arma::uword n = state.gradient_W.n_elem;
  hessian.zeros(n + 1, n + 1);
  hessian(0, 0) = arma::accu(state.curvature_dose);
  const arma::vec cross = state.features.t() * state.curvature_dose;
  hessian.submat(1, 0, n, 0) = cross;
  hessian.submat(0, 1, 0, n) = cross.t();
  hessian.submat(1, 1, n, n) = coefficient;
  return hessian.is_finite();
}

// The same bounded policy proposal is exposed privately for audit fixtures.
// Its feasible slope must be checked after clipping the Newton displacement.
arma::vec direct_feasible_direction(const arma::vec& gradient,
                                    const arma::mat& hessian,
                                    const double alpha,
                                    const DirectPolicy& policy,
                                    std::string& type) {
  arma::vec direction = -gradient;
  direction(0) = policy.project(alpha - gradient(0)) - alpha;
  type = "gradient";
  const bool active = (alpha <= policy.lower && gradient(0) > 0) ||
    (alpha >= policy.upper && gradient(0) < 0);
  const arma::uword first = active ? 1 : 0;
  if (gradient.n_elem <= first || hessian.n_rows != gradient.n_elem ||
      hessian.n_cols != gradient.n_elem || !hessian.is_finite()) return direction;
  arma::mat factor;
  arma::vec intermediate, solution;
  const arma::mat block = hessian.submat(first, first, gradient.n_elem - 1,
                                         gradient.n_elem - 1);
  if (arma::chol(factor, block) &&
      arma::solve(intermediate, arma::trimatl(factor.t()), gradient.tail(gradient.n_elem - first),
                   arma::solve_opts::no_approx) &&
      arma::solve(solution, arma::trimatu(factor), intermediate,
                   arma::solve_opts::no_approx) && solution.is_finite()) {
    arma::vec candidate(gradient.n_elem, arma::fill::zeros);
    candidate.tail(gradient.n_elem - first) = -solution;
    candidate(0) = std::max(policy.lower - alpha,
                            std::min(policy.upper - alpha, candidate(0)));
    const double slope = arma::dot(gradient, candidate);
    if (candidate.is_finite() && std::isfinite(slope) && slope < 0) {
      direction = candidate;
      type = active ? "newton.fixed.intercept" : "newton";
    }
  }
  return direction;
}

arma::vec direct_direction_policy(const DirectEvaluation& state,
                                  const double lambda_reward,
                                  const DirectPolicy& policy,
                                  std::string& type) {
  arma::mat hessian;
  direct_hessian_policy(state, lambda_reward, hessian);
  arma::vec gradient(state.gradient_W.n_elem + 1);
  gradient(0) = state.gradient_alpha;
  gradient.tail(state.gradient_W.n_elem) = state.gradient_W;
  return direct_feasible_direction(gradient, hessian, state.alpha, policy, type);
}

bool direct_cayley(const arma::mat& B, const arma::mat& G,
                    const double step, arma::mat& trial) {
  if (!std::isfinite(step) || step < 0 || !G.is_finite()) return false;
  if (B.n_cols < B.n_rows / 2) {
    const arma::mat U = arma::join_rows(G, B);
    const arma::mat V = arma::join_rows(B, -G);
    const arma::mat system = arma::eye(2 * B.n_cols, 2 * B.n_cols) +
      (0.5 * step) * (V.t() * U);
    arma::mat solution;
    if (!system.is_finite() ||
        !arma::solve(solution, system, V.t() * B, arma::solve_opts::no_approx)) return false;
    trial = B - U * (step * solution);
  } else {
    const arma::mat H = 0.5 * (G * B.t() - B * G.t());
    const arma::mat system = arma::eye(B.n_rows, B.n_rows) + step * H;
    const arma::mat rhs = B - step * (H * B);
    if (!system.is_finite() || !rhs.is_finite() ||
        !arma::solve(trial, system, rhs, arma::solve_opts::no_approx)) return false;
  }
  return trial.is_finite();
}

Rcpp::NumericVector direct_norms(const DirectEvaluation& state) {
  if (state.intercept)
    return Rcpp::NumericVector::create(Rcpp::Named("B") = state.norm_B,
                                       Rcpp::Named("W") = state.norm_W,
                                       Rcpp::Named("intercept") = std::abs(state.residual_alpha),
                                       Rcpp::Named("policy") = state.norm_policy);
  return Rcpp::NumericVector::create(Rcpp::Named("B") = state.norm_B,
                                     Rcpp::Named("W") = state.norm_W);
}

} // namespace

// Private audit access. This uses the production evaluator at arbitrary
// coefficients and does not initialize or optimize a fitted model.
// [[Rcpp::export]]
Rcpp::List pdose_direct_evaluate(const arma::mat B, const arma::mat X,
                                const arma::colvec A, const arma::colvec R,
                                double bw, const arma::colvec W,
                                double lambda_reward, double cayley_step = 0,
                                bool intercept = false,
                                Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds = R_NilValue,
                                double alpha = 0) {
  const DirectPolicy policy = direct_policy(intercept, intercept_bounds, A);
  const DirectEvaluation state = direct_evaluate(B, X, A, R, bw, W, lambda_reward,
                                                true, policy, alpha);
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("valid") = state.valid,
    Rcpp::Named("message") = state.message,
    Rcpp::Named("fn") = state.objective,
    Rcpp::Named("value") = state.value,
    Rcpp::Named("penalty") = state.penalty,
    Rcpp::Named("gB") = state.gradient_B,
    Rcpp::Named("gW") = state.gradient_W,
    Rcpp::Named("gradient.norm") = direct_norms(state),
    Rcpp::Named("kernel") = state.geometry.kernel,
    Rcpp::Named("dose") = state.dose,
    Rcpp::Named("scales") = state.geometry.scales);
  if (state.valid) {
    arma::mat hessian;
    out["hessian.valid"] = direct_hessian_W(state, lambda_reward, hessian);
    out["hessian.W"] = hessian;
    const arma::mat velocity = -(state.gradient_B * B.t() - B * state.gradient_B.t()) * B;
    arma::mat trial;
    out["cayley.valid"] = direct_cayley(B, state.gradient_B, cayley_step, trial);
    out["cayley.B"] = trial;
    out["cayley.velocity"] = velocity;
    out["cayley.slope"] = arma::accu(state.gradient_B % velocity);
    if (intercept) {
      arma::mat policy_hessian;
      out["hessian.policy.valid"] = direct_hessian_policy(state, lambda_reward, policy_hessian);
      out["hessian.policy"] = policy_hessian;
      out["intercept"] = alpha;
      out["intercept.bounds"] = Rcpp::NumericVector::create(policy.lower, policy.upper);
      out["kernel.mean"] = state.geometry.mean.t();
      out["features"] = state.features;
      out["g.intercept"] = state.gradient_alpha;
      out["gradient.intercept"] = state.gradient_alpha;
      std::string type;
      const arma::vec direction = direct_direction_policy(state, lambda_reward, policy, type);
      out["direction.policy"] = direction;
      out["direction.type"] = type;
      out["policy.slope"] = state.gradient_alpha * direction(0) +
        arma::dot(state.gradient_W, direction.tail(W.n_elem));
    }
  }
  return out;
}

// Private regression access to the production projected-Newton proposal.
// [[Rcpp::export]]
Rcpp::List pdose_direct_policy_direction(const arma::vec gradient,
                                        const arma::mat hessian,
                                        double alpha,
                                        const Rcpp::NumericVector intercept_bounds) {
  if (gradient.n_elem < 2 || !gradient.is_finite())
    Rcpp::stop("Policy gradient must be finite with an intercept and coefficients");
  const DirectPolicy policy = direct_policy(true,
    Rcpp::Nullable<Rcpp::NumericVector>(intercept_bounds), arma::vec());
  if (!std::isfinite(alpha) || alpha < policy.lower || alpha > policy.upper)
    Rcpp::stop("Intercept must be finite and within its bounds");
  std::string type;
  const arma::vec direction = direct_feasible_direction(gradient, hessian, alpha, policy, type);
  const double residual = alpha - policy.project(alpha - gradient(0));
  return Rcpp::List::create(Rcpp::Named("direction") = direction,
                            Rcpp::Named("direction.type") = type,
                            Rcpp::Named("slope") = arma::dot(gradient, direction),
                            Rcpp::Named("residual.intercept") = residual);
}

namespace {

struct DirectInitialization {
  arma::vec W;
  arma::vec targets;
  arma::vec gcv;
  double lambda;
  double alpha;
  int count;
};

DirectInitialization direct_initialize(const DirectKernel& geometry,
                                       const arma::vec& A,
                                       const arma::vec& grid,
                                       const arma::vec& R,
                                       const arma::vec& lambda,
                                       const double bw,
                                       const DirectPolicy& policy) {
  DirectInitialization out;
  const arma::uword n = A.n_elem;
  out.count = 1;
  out.alpha = 0;
  out.targets.set_size(n);
  out.gcv.set_size(lambda.n_elem > 1 ? lambda.n_elem : 0);
  const double dose_scale = arma::stddev(A) * bw;
  if (!std::isfinite(dose_scale) || dose_scale <= 0)
    Rcpp::stop("Observed doses must have a positive finite bandwidth scale");
  arma::mat log_dose(n, grid.n_elem);
  for (arma::uword k = 0; k < grid.n_elem; ++k) {
    log_dose.col(k) = -arma::square((A - grid(k)) / dose_scale);
  }
  if (!log_dose.is_finite())
    Rcpp::stop("Initial dose-grid log weights must be finite");
  for (arma::uword i = 0; i < n; ++i) {
    if (i % 32 == 0) Rcpp::checkUserInterrupt();
    double best = -std::numeric_limits<double>::infinity();
    arma::uword selected = 0;
    for (arma::uword k = 0; k < grid.n_elem; ++k) {
      const arma::vec log_weight = geometry.log_kernel.col(i) + log_dose.col(k);
      if (!log_weight.is_finite())
        Rcpp::stop("Initial reward log weights must be finite");
      arma::vec probability = arma::exp(log_weight - log_weight.max());
      probability /= arma::accu(probability);
      const double reward = R(0) + arma::dot(probability, R - R(0));
      if (!std::isfinite(reward)) Rcpp::stop("Initial reward values must be finite");
      if (reward > best) {
        best = reward;
        selected = k;
      }
    }
    out.targets(i) = grid(selected);
  }
  if (policy.intercept) {
    const double target_mean = arma::mean(out.targets);
    out.alpha = policy.project(target_mean);
    out.lambda = lambda(0);
    const bool constant = arma::all(out.targets == out.targets(0));
    if (constant && out.targets(0) >= policy.lower && out.targets(0) <= policy.upper) {
      out.alpha = out.targets(0);
      out.W.zeros(n);
      out.gcv.zeros();
      return out;
    }
    if (lambda.n_elem > 1 && (target_mean <= policy.lower || target_mean >= policy.upper))
      Rcpp::stop("Active intercept.bounds requires a scalar lambda unless the grid targets are an attainable constant");
    // Orthonormal Helmert contrasts remove the intercept null direction.
    // Solving on this space is ordinary intercept KRR, and its trace gives
    // stable residual degrees of freedom without subtracting 1 / lambda.
    arma::mat contrasts(n, n - 1, arma::fill::zeros);
    for (arma::uword j = 0; j + 1 < n; ++j) {
      const double scale = std::sqrt(static_cast<double>((j + 1) * (j + 2)));
      for (arma::uword i = 0; i <= j; ++i) contrasts(i, j) = 1.0 / scale;
      contrasts(j + 1, j) = -static_cast<double>(j + 1) / scale;
    }
    const arma::mat contrast_kernel = arma::symmatu(contrasts.t() * geometry.kernel * contrasts);
    const arma::vec contrast_target = contrasts.t() * (out.targets - target_mean);
    const arma::mat identity = arma::eye(n - 1, n - 1);
    double best_gcv = std::numeric_limits<double>::infinity();
    arma::uword selected = 0;
    if (lambda.n_elem > 1) for (arma::uword k = 0; k < lambda.n_elem; ++k) {
      Rcpp::checkUserInterrupt();
      arma::mat inverse;
      if (!arma::solve(inverse, contrast_kernel + lambda(k) * identity, identity,
                       arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
          !inverse.is_finite())
        Rcpp::stop("Initial intercept kernel ridge system could not be solved");
      const double denominator = arma::trace(inverse);
      const double ratio = arma::norm(inverse * contrast_target, 2) / denominator;
      out.gcv(k) = static_cast<double>(n) * ratio * ratio;
      if (!std::isfinite(out.gcv(k)) || denominator <= 0)
        Rcpp::stop("Initial intercept kernel ridge GCV score must be finite");
      if (out.gcv(k) < best_gcv) {
        best_gcv = out.gcv(k);
        selected = k;
      }
    }
    out.lambda = lambda(selected);
    arma::vec solution;
    if (!arma::solve(solution, contrast_kernel + out.lambda * identity, contrast_target,
                     arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
        !solution.is_finite())
      Rcpp::stop("Initial intercept kernel ridge coefficients could not be solved");
    out.W = contrasts * solution;
    return out;
  }
  const arma::mat identity = arma::eye(n, n);
  if (lambda.n_elem == 1) {
    out.lambda = lambda(0);
    if (!arma::solve(out.W, geometry.kernel + out.lambda * identity, out.targets,
                     arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
        !out.W.is_finite())
      Rcpp::stop("Initial kernel ridge coefficients could not be solved");
    return out;
  }
  double best_gcv = std::numeric_limits<double>::infinity();
  arma::uword selected = 0;
  for (arma::uword k = 0; k < lambda.n_elem; ++k) {
    Rcpp::checkUserInterrupt();
    const arma::mat system = geometry.kernel + lambda(k) * identity;
    arma::mat solution;
    if (!arma::solve(solution, system, identity,
                     arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
        !solution.is_finite())
      Rcpp::stop("Initial kernel ridge system could not be solved");
    // I - K (K + lambda I)^(-1) = lambda (K + lambda I)^(-1).
    // Cancel lambda from numerator and denominator before squaring. This
    // is the same GCV criterion without subtracting nearly equal matrices.
    const double denominator = arma::trace(solution);
    const double ratio = arma::norm(solution * out.targets, 2) / denominator;
    out.gcv(k) = static_cast<double>(n) * ratio * ratio;
    if (!std::isfinite(out.gcv(k)) || denominator <= 0)
      Rcpp::stop("Initial kernel ridge GCV score must be finite");
    // Strict comparison retains the first candidate in an exact tie.
    if (out.gcv(k) < best_gcv) {
      best_gcv = out.gcv(k);
      selected = k;
    }
  }
  out.lambda = lambda(selected);
  if (!arma::solve(out.W, geometry.kernel + out.lambda * identity, out.targets,
                   arma::solve_opts::likely_sympd + arma::solve_opts::no_approx) ||
      !out.W.is_finite())
    Rcpp::stop("Initial kernel ridge coefficients could not be solved");
  return out;
}

struct DirectControls {
  double rho, eta, tau, btol, ftol, gtol;
  int maxitr, verbose, B_trial_budget;
};

void direct_validate(const arma::mat& B, const arma::mat& X,
                     const arma::vec& A, const arma::vec& R,
                     const arma::mat& a_dist, const arma::vec& grid,
                     const arma::vec& lambda, const double bw,
                     const double lambda_reward, const DirectControls& control) {
  if (A.n_elem != X.n_rows || R.n_elem != X.n_rows ||
      !A.is_finite() || !R.is_finite())
    Rcpp::stop("Direct-dose observations must be finite with matching lengths");
  if (grid.n_elem == 0 || !grid.is_finite() || a_dist.n_rows != X.n_rows ||
      a_dist.n_cols != grid.n_elem)
    Rcpp::stop("Initial dose grid has invalid dimensions or values");
  if (lambda.n_elem == 0 || !lambda.is_finite() || arma::any(lambda <= 0))
    Rcpp::stop("lambda must contain positive finite initialization ridge values");
  if (!std::isfinite(lambda_reward) || lambda_reward < 0)
    Rcpp::stop("lambda.reward must be a finite nonnegative scalar");
  if (!std::isfinite(control.rho) || control.rho <= 0 || control.rho >= 1 ||
      !std::isfinite(control.eta) || control.eta <= 0 || control.eta >= 1 ||
      !std::isfinite(control.tau) || control.tau <= 0 ||
      !std::isfinite(control.btol) || control.btol <= 0 ||
      !std::isfinite(control.ftol) || control.ftol <= 0 ||
      !std::isfinite(control.gtol) || control.gtol <= 0 || control.maxitr < 0 ||
      control.B_trial_budget < 0 || control.B_trial_budget > 30)
    Rcpp::stop("Direct-dose optimizer controls are invalid");
  DirectKernel geometry;
  std::string message;
  if (!direct_kernel(B, X, bw, geometry, message)) Rcpp::stop(message);
  if (arma::norm(B.t() * B - arma::eye(B.n_cols, B.n_cols), "fro") > 1e-8)
    Rcpp::stop("Initial direct-dose basis must have orthonormal columns");
}

void direct_trace(Rcpp::List& trace, const bool collect,
                  const std::string& event, const std::string& block,
                  const int sweep, const arma::mat& B_start,
                  const arma::vec& W_start, const arma::mat& B,
                  const arma::vec& W, const double fn_start,
                  const DirectEvaluation& state, const double slope,
                  const double step, const bool accepted,
                  const std::string& reason = "",
                  const std::string& direction_type = "",
                  const arma::vec& direction_W = arma::vec(),
                  const double alpha_start = NA_REAL,
                  const double direction_alpha = NA_REAL) {
  if (!collect) return;
  Rcpp::List entry = Rcpp::List::create(
    Rcpp::Named("event") = event, Rcpp::Named("block") = block,
    Rcpp::Named("sweep") = sweep,
    Rcpp::Named("B.start") = B_start, Rcpp::Named("W.start") = W_start,
    Rcpp::Named("B") = B, Rcpp::Named("W") = W,
    Rcpp::Named("fn.start") = fn_start, Rcpp::Named("fn") = state.objective,
    Rcpp::Named("value") = state.value, Rcpp::Named("penalty") = state.penalty,
    Rcpp::Named("slope") = slope, Rcpp::Named("step") = step,
    Rcpp::Named("accepted") = accepted, Rcpp::Named("finite") = state.valid,
    Rcpp::Named("reason") = reason, Rcpp::Named("message") = state.message,
    Rcpp::Named("direction.type") = direction_type,
    Rcpp::Named("direction.W") = direction_W);
  if (state.intercept) {
    entry["intercept.start"] = alpha_start;
    entry["intercept"] = state.alpha;
    entry["direction.intercept"] = direction_alpha;
    entry["kernel.mean"] = state.geometry.mean.t();
    entry["gradient.intercept"] = state.gradient_alpha;
    entry["gradient.norm"] = direct_norms(state);
  }
  trace.push_back(entry);
}

double direct_next_step(const double accepted, const double eta,
                         const double cap = 1e10) {
  return accepted >= cap * eta ? cap : accepted / eta;
}

Rcpp::List direct_fit(arma::mat B, const arma::mat& X, const arma::vec& A,
                      const arma::mat& a_dist, const arma::vec& grid,
                      const arma::vec& R, const arma::vec& lambda,
                      const double bw, const double lambda_reward,
                      const DirectControls& control, const bool collect,
                      const Rcpp::Nullable<Rcpp::NumericVector>& W_initial,
                      const DirectPolicy& policy,
                      const Rcpp::Nullable<Rcpp::NumericVector>& alpha_initial) {
  direct_validate(B, X, A, R, a_dist, grid, lambda, bw, lambda_reward, control);
  if (policy.intercept && lambda_reward <= 0)
    Rcpp::stop("lambda.reward must be positive when intercept = TRUE");
  if (!policy.intercept && alpha_initial.isNotNull())
    Rcpp::stop("alpha_initial requires intercept = TRUE");
  DirectKernel initial_geometry;
  std::string message;
  if (!direct_kernel(B, X, bw, initial_geometry, message)) Rcpp::stop(message);
  DirectInitialization initialization;
  if (W_initial.isNotNull()) {
    initialization.W = Rcpp::as<arma::vec>(W_initial.get());
    initialization.lambda = lambda(0);
    initialization.count = 0;
    initialization.alpha = policy.intercept ? policy.project(arma::mean(A)) : 0;
  } else {
    initialization = direct_initialize(initial_geometry, A, grid, R, lambda, bw, policy);
  }
  if (alpha_initial.isNotNull()) {
    const Rcpp::NumericVector initial(alpha_initial.get());
    if (initial.size() != 1 || !std::isfinite(initial[0]) ||
        initial[0] < policy.lower || initial[0] > policy.upper)
      Rcpp::stop("alpha_initial must be a finite scalar within intercept.bounds");
    initialization.alpha = initial[0];
  }
  arma::vec W = initialization.W;
  double alpha = initialization.alpha;
  DirectEvaluation state = direct_evaluate(B, X, A, R, bw, W, lambda_reward, true, policy, alpha);
  if (!state.valid) Rcpp::stop("Initial direct-dose evaluation failed: " + state.message);
  Rcpp::List trace;
  if (collect) {
    direct_trace(trace, true, "initialize", "", 0, B, W, B, W,
                  state.objective, state, NA_REAL, 0, true, "", "", arma::vec(), alpha);
    Rcpp::List entry = trace[0];
    entry["initialization.count"] = initialization.count;
    entry["lambda.init"] = initialization.lambda;
    entry["grid.targets"] = initialization.targets;
    entry["GCV"] = initialization.gcv;
    trace[0] = entry;
  }
  int completed = 0;
  int small_sweeps = 0;
  double W_step = std::min(control.tau, 1.0);
  double B_step = std::min(control.tau, 1e10);
  std::string reason = "iteration_limit";
  if (control.verbose > 1)
    Rcpp::Rcout << "Initial direct-dose objective: " << state.objective << std::endl;

  while (completed < control.maxitr) {
    Rcpp::checkUserInterrupt();
    if (state.norm_B <= control.gtol && state.norm_policy <= control.gtol) {
      reason = "gradient_tolerance";
      break;
    }
    const arma::mat previous_B = B;
    const arma::vec previous_W = W;
    const double previous_alpha = alpha;
    const double previous_objective = state.objective;
    arma::mat proposed_B = B;
    arma::vec proposed_W = W;
    double proposed_alpha = alpha;
    DirectEvaluation proposed = state;
    bool failed = false;
    std::string failed_block;

    const std::string policy_block = policy.intercept ? "policy" : "W";
    if (proposed.norm_policy > control.gtol) {
      failed_block = policy_block;
      const arma::vec gradient = proposed.gradient_W;
      std::string direction_type;
      arma::vec direction;
      double direction_alpha = 0;
      if (policy.intercept) {
        const arma::vec joint = direct_direction_policy(proposed, lambda_reward, policy, direction_type);
        direction_alpha = joint(0);
        direction = joint.tail(W.n_elem);
      } else direction = direct_direction_W(proposed, lambda_reward, direction_type);
      const double slope = arma::dot(gradient, direction) +
        (policy.intercept ? proposed.gradient_alpha * direction_alpha : 0);
      if (!std::isfinite(slope) || slope >= 0) {
        reason = std::isfinite(slope) ? "stagnation" : "numerical_failure";
        failed = true;
      } else {
        bool accepted = false;
        double step = W_step;
        for (int trial = 0; trial < 30; ++trial) {
          Rcpp::checkUserInterrupt();
          const arma::vec trial_W = W + step * direction;
          // The displacement is feasible and step <= 1. Projection here
          // only absorbs rounding when the full displacement hits a bound.
          const double trial_alpha = policy.intercept ?
            policy.project(alpha + step * direction_alpha) : 0;
          DirectEvaluation candidate;
          if (trial_W.is_finite())
            candidate = direct_evaluate(B, X, A, R, bw, trial_W, lambda_reward,
                                         false, policy, trial_alpha);
          else {
            candidate.intercept = policy.intercept;
            candidate.alpha = trial_alpha;
            candidate.message = "Coefficient trial is nonfinite";
          }
          const bool moved = trial_W.n_elem == W.n_elem &&
            (arma::any(trial_W != W) || trial_alpha != alpha);
          const double threshold = proposed.objective + control.rho * step * slope;
          accepted = moved && candidate.valid && std::isfinite(threshold) &&
            candidate.objective <= threshold;
          direct_trace(trace, collect, "trial", policy_block, completed + 1, B, W,
                        B, trial_W, proposed.objective, candidate, slope, step, accepted,
                        "", direction_type, direction, alpha, direction_alpha);
          if (candidate.valid && !moved) {
            reason = "stagnation";
            failed = true;
            break;
          }
          if (accepted) {
            proposed_W = trial_W;
            proposed_alpha = trial_alpha;
            proposed = candidate;
            W_step = direct_next_step(step, control.eta, 1.0);
            break;
          }
          step *= control.eta;
        }
        if (!accepted && !failed) {
          reason = policy.intercept ? "policy_line_search_failed" : "w_line_search_failed";
          failed = true;
        }
      }
    } else {
      direct_trace(trace, collect, "skip", policy_block, completed + 1, B, W, B, W,
                    proposed.objective, proposed, 0, 0, true, "", "", arma::vec(), alpha);
    }

    // Recompute the direction derivative at provisional coefficients. A
    // derivative failure rolls back their accepted objective-only trial.
    if (!failed) {
      proposed = direct_evaluate(proposed_B, X, A, R, bw, proposed_W, lambda_reward,
                                   true, policy, proposed_alpha);
      if (!proposed.valid) {
        reason = "numerical_failure";
        failed_block = "B";
        failed = true;
      }
    }
    if (!failed && proposed.norm_B > control.gtol) {
      failed_block = "B";
      const arma::mat gradient = proposed.gradient_B;
      const arma::mat velocity = -(gradient * B.t() - B * gradient.t()) * B;
      const double slope = arma::accu(gradient % velocity);
      if (!std::isfinite(slope) || slope >= 0) {
        reason = std::isfinite(slope) ? "stagnation" : "numerical_failure";
        failed = true;
      } else {
        bool accepted = false;
        double step = B_step;
        for (int trial = 0; trial < control.B_trial_budget; ++trial) {
          Rcpp::checkUserInterrupt();
          arma::mat trial_B;
          DirectEvaluation candidate;
          const bool feasible = direct_cayley(B, gradient, step, trial_B);
          if (feasible)
            candidate = direct_evaluate(trial_B, X, A, R, bw, proposed_W, lambda_reward,
                                         false, policy, proposed_alpha);
          else {
            candidate.intercept = policy.intercept;
            candidate.alpha = proposed_alpha;
            candidate.message = "Cayley trial is nonfinite or unsolved";
          }
          const bool moved = feasible && arma::accu(trial_B != B) > 0;
          const double threshold = proposed.objective + control.rho * step * slope;
          accepted = moved && candidate.valid && std::isfinite(threshold) &&
            candidate.objective <= threshold;
          direct_trace(trace, collect, "trial", "B", completed + 1, B, proposed_W,
                        trial_B, proposed_W, proposed.objective, candidate, slope, step, accepted,
                        "", "", arma::vec(), proposed_alpha);
          if (candidate.valid && !moved) {
            reason = "stagnation";
            failed = true;
            break;
          }
          if (accepted) {
            proposed_B = trial_B;
            proposed = candidate;
            B_step = direct_next_step(step, control.eta);
            break;
          }
          step *= control.eta;
        }
        if (!accepted && !failed) {
          reason = "B_line_search_failed";
          failed = true;
        }
      }
    } else if (!failed) {
      direct_trace(trace, collect, "skip", "B", completed + 1, B, proposed_W,
                    B, proposed_W, proposed.objective, proposed, 0, 0, true,
                    "", "", arma::vec(), proposed_alpha);
    }

    // Both fresh gradients and all components belong to the final pair,
    // including at the last budgeted sweep, before the atomic commit.
    if (!failed) {
      proposed = direct_evaluate(proposed_B, X, A, R, bw, proposed_W, lambda_reward,
                                   true, policy, proposed_alpha);
      if (!proposed.valid) {
        reason = "numerical_failure";
        failed_block = "B";
        failed = true;
      }
    }
    if (failed) {
      direct_trace(trace, collect, "rollback", failed_block, completed + 1,
                    proposed_B, proposed_W, B, W, previous_objective, state,
                    NA_REAL, 0, false, reason, "", arma::vec(), proposed_alpha);
      break;
    }
    B = proposed_B;
    W = proposed_W;
    alpha = proposed_alpha;
    state = proposed;
    ++completed;
    direct_trace(trace, collect, "commit", "", completed, previous_B,
                  previous_W, B, W, previous_objective, state, NA_REAL, 0, true,
                  "", "", arma::vec(), previous_alpha);
    if (state.norm_B <= control.gtol && state.norm_policy <= control.gtol) {
      reason = "gradient_tolerance";
      break;
    }
    const double B_change = arma::norm(B - previous_B, "fro") /
      std::max(1.0, arma::norm(previous_B, "fro"));
    const double W_change = arma::norm(W - previous_W, 2) /
      std::max(1.0, arma::norm(previous_W, 2));
    const double alpha_change = std::abs(alpha - previous_alpha) /
      std::max(1.0, std::abs(previous_alpha));
    const double objective_change = std::abs(state.objective - previous_objective) /
      (1.0 + std::abs(previous_objective));
    if (B_change <= control.btol && W_change <= control.btol &&
        (!policy.intercept || alpha_change <= control.btol) &&
        objective_change <= control.ftol) ++small_sweeps;
    else small_sweeps = 0;
    if (small_sweeps >= 5) {
      reason = "stagnation";
      break;
    }
    if (control.verbose > 1 && completed % 10 == 0)
      Rcpp::Rcout << "At sweep " << completed << ", objective = " << state.objective << std::endl;
  }
  // Also covers an initially stationary pair with a zero sweep budget.
  if (state.norm_B <= control.gtol && state.norm_policy <= control.gtol)
    reason = "gradient_tolerance";
  const DirectEvaluation fresh = direct_evaluate(B, X, A, R, bw, W, lambda_reward,
                                                true, policy, alpha);
  if (fresh.valid) state = fresh;
  else {
    reason = "numerical_failure";
    state.norm_B = NA_REAL;
    state.norm_W = NA_REAL;
    state.residual_alpha = NA_REAL;
    state.norm_policy = NA_REAL;
  }
  const bool converged = reason == "gradient_tolerance";
  if (control.verbose > 0) {
    Rcpp::Rcout << "Direct-dose stop: " << reason << "; committed sweeps: " << completed <<
      "; objective: " << state.objective << "; gradient norms (B, W): " <<
      state.norm_B << ", " << state.norm_W << std::endl;
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("B") = B, Rcpp::Named("W") = W,
    Rcpp::Named("fn") = state.objective, Rcpp::Named("itr") = completed,
    Rcpp::Named("converge") = converged,
    Rcpp::Named("value") = state.value, Rcpp::Named("penalty") = state.penalty,
    Rcpp::Named("lambda.init") = initialization.lambda,
    Rcpp::Named("lambda.reward") = lambda_reward,
    Rcpp::Named("gradient.norm") = direct_norms(state),
    Rcpp::Named("stop.reason") = reason,
    Rcpp::Named("objective.version") = policy.intercept ? "direct_joint_intercept_v2" : "direct_joint_v1");
  if (policy.intercept) {
    out["intercept"] = alpha;
    out["intercept.bounds"] = Rcpp::NumericVector::create(policy.lower, policy.upper);
    out["kernel.mean"] = state.geometry.mean.t();
    out["gradient.intercept"] = state.gradient_alpha;
  }
  if (collect) out["trace"] = trace;
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List pdose_direct_solver(arma::mat B,
                              const arma::mat X,
                              const arma::colvec A,
                              const arma::mat a_dist,
                              const arma::colvec a_seq,
                              const arma::colvec R,
                              const arma::colvec lambda,
                              double bw, double rho, double eta, double gamma,
                              double tau, double epsilon, double btol,
                              double ftol, double gtol, int maxitr,
                              int verbose, int ncore, double lambda_reward = 0.1,
                              bool intercept = false,
                              Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds = R_NilValue) {
  // gamma and epsilon are retained only for internal calling compatibility.
  (void)gamma;
  (void)epsilon;
  (void)ncore;
  const DirectControls control = {rho, eta, tau, btol, ftol, gtol,
                                  maxitr, verbose, 30};
  return direct_fit(B, X, A, a_dist, a_seq, R, lambda, bw, lambda_reward,
                     control, false, R_NilValue,
                     direct_policy(intercept, intercept_bounds, A), R_NilValue);
}

// Private audit access to the same loop. A zero direction trial budget
// exercises atomic rollback; production always uses thirty trials.
// [[Rcpp::export]]
Rcpp::List pdose_direct_audit_solver(arma::mat B,
                                    const arma::mat X,
                                    const arma::colvec A,
                                    const arma::mat a_dist,
                                    const arma::colvec a_seq,
                                    const arma::colvec R,
                                    const arma::colvec lambda,
                                    double bw, double rho, double eta, double gamma,
                                    double tau, double epsilon, double btol,
                                    double ftol, double gtol, int maxitr,
                                    int verbose, int ncore, double lambda_reward = 0.1,
                                    bool trace = true, int B_trial_budget = 30,
                                    Rcpp::Nullable<Rcpp::NumericVector> W_initial = R_NilValue,
                                    bool intercept = false,
                                    Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds = R_NilValue,
                                    Rcpp::Nullable<Rcpp::NumericVector> alpha_initial = R_NilValue) {
  (void)gamma;
  (void)epsilon;
  (void)ncore;
  const DirectControls control = {rho, eta, tau, btol, ftol, gtol,
                                  maxitr, verbose, B_trial_budget};
  return direct_fit(B, X, A, a_dist, a_seq, R, lambda, bw, lambda_reward,
                     control, trace, W_initial,
                     direct_policy(intercept, intercept_bounds, A), alpha_initial);
}
