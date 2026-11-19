// Copyright orthoDr authors. Distributed under GPL (>= 2).
#include <RcppArmadillo.h>
#include "utilities.h"
#include "orthoDr_pdose.h"
#include "solver_numerics.h"
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

namespace {

const double pseudo_feasibility_tolerance = 1e-8;
const double pseudo_final_dose_factor = 0.70710678118654752440;

// Center before squaring and scale the sum of squares. Long-double sums also
// avoid manufacturing variation when many identical decimal values are added.
bool pseudo_center(const arma::vec& x, arma::vec& centered, double& mean,
                   double& scale) {
  if (x.n_elem < 2 || !x.is_finite()) return false;
  long double offset = 0;
  const long double anchor = x(0);
  for (arma::uword i = 0; i < x.n_elem; ++i)
    offset += static_cast<long double>(x(i)) - anchor;
  const long double location = anchor + offset / x.n_elem;
  long double largest = 0;
  for (arma::uword i = 0; i < x.n_elem; ++i)
    largest = std::max(largest, std::abs(static_cast<long double>(x(i)) - location));
  if (!std::isfinite(location) || !std::isfinite(largest) || largest == 0) return false;
  long double squares = 0;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    const long double z = (static_cast<long double>(x(i)) - location) / largest;
    squares += z * z;
  }
  const long double deviation = largest * std::sqrt(squares / (x.n_elem - 1));
  mean = static_cast<double>(location);
  scale = static_cast<double>(deviation);
  if (!std::isfinite(mean) || !std::isfinite(scale) || scale <= 0) return false;
  centered.set_size(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i)
    centered(i) = static_cast<double>(static_cast<long double>(x(i)) - location);
  return centered.is_finite();
}

bool pseudo_constant(const arma::vec& x) {
  for (arma::uword i = 1; i < x.n_elem; ++i)
    if (x(i) != x(0)) return false;
  return x.n_elem > 0;
}

struct PseudoData {
  arma::mat working_X;
  arma::vec q, A;
  double reward_scale, dose_scale, bw;
};

PseudoData pseudo_data(const arma::mat& X, const arma::vec& R,
                       const arma::vec& A, const double bw) {
  if (X.n_rows < 2 || X.n_cols == 0 || R.n_elem != X.n_rows ||
      A.n_elem != X.n_rows || !X.is_finite() || !R.is_finite() ||
      !A.is_finite() || !std::isfinite(bw) || bw <= 0)
    Rcpp::stop("Pseudo-direct data and bandwidth must be finite with valid dimensions and a positive bandwidth");
  if (pseudo_constant(R))
    Rcpp::stop("Constant rewards do not identify directions or a preferred dose");
  PseudoData out;
  out.bw = bw;
  out.A = A;
  arma::vec centered;
  double location, scale;
  if (!pseudo_center(R, centered, location, out.reward_scale))
    Rcpp::stop("Varying rewards must have a positive finite representable sample standard deviation");
  out.q = centered / out.reward_scale;
  if (!out.q.is_finite()) Rcpp::stop("Normalized rewards must be finite");
  if (!pseudo_center(A, centered, location, scale))
    Rcpp::stop("Observed doses must have a positive finite sample standard deviation");
  out.dose_scale = scale * bw;
  if (!std::isfinite(out.dose_scale) || out.dose_scale <= 0)
    Rcpp::stop("Observed dose bandwidth scale must be positive and finite");
  // Center the projected index below. Subtracting separate feature means
  // here can discard a small contrast even if the index is well resolved.
  out.working_X = X.each_row() - dose_origin(X);
  return out;
}

struct PseudoEvaluation {
  bool valid;
  std::string message;
  double objective, norm;
  arma::mat gradient, residual, log_kernel, kernel;
  arma::vec fitted;
  PseudoEvaluation() : valid(false), objective(NA_REAL), norm(NA_REAL) {}
};

PseudoEvaluation pseudo_evaluate(const arma::mat& B, const PseudoData& data,
                                 const bool derivatives) {
  PseudoEvaluation out;
  const arma::uword n = data.A.n_elem;
  if (B.n_rows != data.working_X.n_cols || B.n_cols == 0 ||
      B.n_cols > B.n_rows || !B.is_finite()) {
    out.message = "Pseudo-direct basis dimensions and entries must be valid and finite";
    return out;
  }
  arma::mat projected = data.working_X * B;
  arma::mat centered(n, B.n_cols);
  arma::rowvec scales(B.n_cols);
  out.log_kernel.zeros(n, n);
  for (arma::uword k = 0; k < B.n_cols; ++k) {
    Rcpp::checkUserInterrupt();
    arma::vec column;
    double location;
    if (!pseudo_center(projected.col(k), column, location, scales(k))) {
      out.message = "Projected covariates must have positive finite sample standard deviations";
      return out;
    }
    centered.col(k) = column;
    const double scale = scales(k) * data.bw;
    if (!std::isfinite(scale) || scale <= 0) {
      out.message = "Projected bandwidth scales must be positive and finite";
      return out;
    }
    const arma::vec z = column / scale;
    const arma::mat delta = arma::repmat(z, 1, n) - arma::repmat(z.t(), n, 1);
    out.log_kernel -= 0.5 * arma::square(delta);
  }
  if (!out.log_kernel.is_finite()) {
    out.message = "Covariate log-kernel arithmetic must be finite";
    return out;
  }
  out.kernel = arma::exp(out.log_kernel);
  arma::mat sensitivity;
  if (derivatives) sensitivity.set_size(n, n);
  out.fitted.set_size(n);
  out.objective = 0;
  for (arma::uword i = 0; i < n; ++i) {
    if (i % 32 == 0) Rcpp::checkUserInterrupt();
    const arma::vec delta = (data.A - data.A(i)) / data.dose_scale;
    const arma::vec log_weight = out.log_kernel.col(i) - 0.5 * arma::square(delta);
    if (!log_weight.is_finite()) {
      out.message = "Training joint log-weight arithmetic must be finite";
      return out;
    }
    arma::vec probability = arma::exp(log_weight - log_weight.max());
    const double denominator = arma::accu(probability);
    if (!std::isfinite(denominator) || denominator <= 0) {
      out.message = "Training reward weights cannot be normalized";
      return out;
    }
    probability /= denominator;
    out.fitted(i) = arma::dot(probability, data.q);
    const double residual = out.fitted(i) - data.q(i);
    out.objective += residual * residual / static_cast<double>(n);
    if (derivatives)
      sensitivity.col(i) = (2.0 * residual / static_cast<double>(n)) *
        probability % (data.q - out.fitted(i));
  }
  if (!std::isfinite(out.objective)) {
    out.message = "Normalized pseudo-direct objective must be finite";
    return out;
  }
  if (derivatives) {
    const arma::mat symmetric = sensitivity + sensitivity.t();
    arma::mat gradient_projected(n, B.n_cols);
    for (arma::uword k = 0; k < B.n_cols; ++k) {
      const double scale = scales(k);
      const arma::vec z = centered.col(k) / (data.bw * scale);
      const arma::mat delta = arma::repmat(z, 1, n) - arma::repmat(z.t(), n, 1);
      // The second term differentiates b_k' cov(X) b_k, with divisor n-1.
      const double scale_sensitivity = arma::accu(sensitivity % arma::square(delta));
      gradient_projected.col(k) =
        -arma::sum(symmetric % delta, 1) / (data.bw * scale) +
        scale_sensitivity * (centered.col(k) / scale) /
          ((static_cast<double>(n) - 1.0) * scale);
    }
    out.gradient = data.working_X.t() * gradient_projected;
    out.residual = out.gradient - B * (out.gradient.t() * B);
    out.norm = arma::norm(out.residual, "fro");
    if (!out.gradient.is_finite() || !out.residual.is_finite() || !std::isfinite(out.norm)) {
      out.message = "Normalized pseudo-direct derivatives must be finite";
      return out;
    }
  }
  out.valid = true;
  return out;
}

double pseudo_feasibility(const arma::mat& B) {
  return arma::norm(B.t() * B - arma::eye(B.n_cols, B.n_cols), "fro");
}

struct PseudoCurve {
  bool dense;
  arma::mat H, HB, U, VU, VB;
  double slope;
  PseudoCurve(const arma::mat& B, const arma::mat& G, const bool use_dense)
    : dense(use_dense) {
    if (dense) {
      const arma::mat GB = G * B.t();
      H = 0.5 * (GB - GB.t());
      HB = H * B;
      slope = 2.0 * arma::accu(G % HB);
    } else {
      U = arma::join_rows(G, B);
      const arma::mat V = arma::join_rows(B, -G);
      VU = V.t() * U;
      VB = V.t() * B;
      slope = arma::accu(G % (U * VB));
    }
  }
  bool trial(const arma::mat& B, const double step, arma::mat& result) const {
    bool solved;
    if (dense) {
      const arma::mat lhs = arma::eye(B.n_rows, B.n_rows) + step * H;
      const arma::mat rhs = B - step * HB;
      solved = arma::solve(result, lhs, rhs, arma::solve_opts::no_approx);
    } else {
      arma::mat coefficients;
      const arma::mat lhs = arma::eye(VU.n_rows, VU.n_cols) + 0.5 * step * VU;
      solved = arma::solve(coefficients, lhs, VB, arma::solve_opts::no_approx);
      if (solved) result = B - U * (step * coefficients);
    }
    return solved && result.is_finite() &&
      pseudo_feasibility(result) <= pseudo_feasibility_tolerance;
  }
};

struct PseudoTargets {
  arma::vec dose;
  arma::mat scores;
};

struct CompensatedSum {
  double sum = 0;
  double correction = 0;

  void add(const double value) {
    const double next = sum + value;
    // Recover the rounding error even when the new summand is larger.
    if (std::abs(sum) >= std::abs(value))
      correction += (sum - next) + value;
    else
      correction += (value - next) + sum;
    sum = next;
  }

  double total() const { return sum + correction; }
};

PseudoTargets pseudo_targets(const PseudoEvaluation& state, const PseudoData& data,
                              const arma::vec& grid, const double dose_factor,
                              const bool keep_scores) {
  if (grid.n_elem == 0 || !grid.is_finite() ||
      !std::isfinite(dose_factor) || dose_factor <= 0)
    Rcpp::stop("Dose grid and final dose bandwidth factor must be finite and valid");
  const double scale = data.dose_scale * dose_factor;
  if (!std::isfinite(scale) || scale <= 0)
    Rcpp::stop("Final dose bandwidth scale must be positive and finite");
  const arma::uword n = data.A.n_elem;
  PseudoTargets out;
  out.dose.set_size(n);
  if (keep_scores) out.scores.set_size(n, grid.n_elem);
  arma::mat log_dose(n, grid.n_elem);
  for (arma::uword j = 0; j < grid.n_elem; ++j) {
    // Form dose contributions from raw doses. Exponentiated a_dist loses the
    // information needed when every product in a distant grid cell underflows.
    const arma::vec delta = (data.A - grid(j)) / scale;
    log_dose.col(j) = -0.5 * arma::square(delta);
  }
  if (!log_dose.is_finite())
    Rcpp::stop("Dose extraction has invalid joint log-weight geometry");
  const double radius_factor = 4.0 * std::numeric_limits<double>::epsilon();
  arma::vec upper(grid.n_elem);
  for (arma::uword i = 0; i < n; ++i) {
    if (i % 32 == 0) Rcpp::checkUserInterrupt();
    double best_lower = -std::numeric_limits<double>::infinity();
    for (arma::uword j = 0; j < grid.n_elem; ++j) {
      const arma::vec log_weight = state.log_kernel.col(i) + log_dose.col(j);
      if (!log_weight.is_finite())
        Rcpp::stop("Dose extraction has no finite joint log-weight reference");
      const arma::vec weights = arma::exp(log_weight - log_weight.max());
      CompensatedSum weight_sum, reward_sum, magnitude_sum;
      for (arma::uword k = 0; k < n; ++k) {
        const double product = weights(k) * data.q(k);
        weight_sum.add(weights(k));
        reward_sum.add(product);
        magnitude_sum.add(std::abs(product));
      }
      const double denominator = weight_sum.total();
      if (!std::isfinite(denominator) || denominator <= 0)
        Rcpp::stop("Dose extraction weights cannot be normalized");
      const double score = reward_sum.total() / denominator;
      if (!std::isfinite(score)) Rcpp::stop("Dose-grid reward scores must be finite");
      if (keep_scores) out.scores(i, j) = score;
      const double radius = radius_factor * (magnitude_sum.total() / denominator);
      best_lower = std::max(best_lower, score - radius);
      upper(j) = score + radius;
    }
    // At least the candidate attaining best_lower reaches that endpoint.
    arma::uword selected = 0;
    while (upper(selected) < best_lower) ++selected;
    out.dose(i) = grid(selected);
  }
  return out;
}

bool pseudo_accurate_solve(const arma::mat& C, const arma::mat& rhs,
                           const arma::mat& solution) {
  if (!solution.is_finite()) return false;
  const double residual = arma::norm(C * solution - rhs, "inf");
  const long double denominator = static_cast<long double>(arma::norm(C, "inf")) *
    arma::norm(solution, "inf") + arma::norm(rhs, "inf");
  return std::isfinite(residual) && std::isfinite(denominator) &&
    (denominator == 0 ? residual == 0 :
      residual / denominator <= 128.0 * C.n_rows * std::numeric_limits<double>::epsilon());
}

bool pseudo_chol_solve(const arma::mat& C, const arma::mat& factor,
                       const arma::mat& rhs, arma::mat& solution) {
  arma::mat intermediate;
  return arma::solve(intermediate, arma::trimatl(factor.t()), rhs,
                      arma::solve_opts::no_approx) &&
    arma::solve(solution, arma::trimatu(factor), intermediate,
                  arma::solve_opts::no_approx) &&
    pseudo_accurate_solve(C, rhs, solution);
}

struct PseudoRidge {
  arma::vec W, scores;
  arma::uvec valid;
  double alpha, lambda;
};

PseudoRidge pseudo_ridge(const arma::mat& K, const arma::vec& target,
                         const arma::vec& lambda, const bool intercept) {
  if (lambda.n_elem == 0 || !lambda.is_finite() || arma::any(lambda <= 0))
    Rcpp::stop("Pseudo-direct lambda must contain positive finite ridge values");
  const arma::uword n = target.n_elem;
  const arma::mat identity = arma::eye(n, n);
  const arma::mat ones(n, 1, arma::fill::ones);
  const bool constant = pseudo_constant(target);
  double target_mean = 0, unused_scale;
  arma::vec centered_target;
  if (intercept && !constant && !pseudo_center(target, centered_target, target_mean, unused_scale))
    Rcpp::stop("Pseudo-direct ridge targets cannot be centered with finite arithmetic");
  PseudoRidge out;
  out.alpha = 0;
  out.lambda = NA_REAL;
  out.scores.set_size(lambda.n_elem);
  out.scores.fill(NA_REAL);
  out.valid.zeros(lambda.n_elem);
  double best = std::numeric_limits<double>::infinity();
  std::vector<double> skipped;
  for (arma::uword j = 0; j < lambda.n_elem; ++j) {
    Rcpp::checkUserInterrupt();
    arma::mat C = K;
    C.diag() += lambda(j);
    arma::mat factor, inverse, v, solution;
    bool valid = C.is_finite() && arma::chol(factor, C) &&
      pseudo_chol_solve(C, factor, identity, inverse);
    arma::vec W;
    double alpha = 0;
    long double denominator = valid ? arma::trace(inverse) : 0;
    if (valid && intercept) {
      valid = pseudo_chol_solve(C, factor, ones, v);
      const double t = valid ? arma::accu(v) : 0;
      valid = valid && std::isfinite(t) && t > 0;
      if (valid) {
        // Divide before multiplying: squaring v first underflows at large
        // finite ridge values, even when its degrees-of-freedom term is valid.
        const long double vnorm = arma::norm(v, 2);
        denominator -= vnorm * (vnorm / t);
        if (constant) {
          alpha = target(0);
          W.zeros(n);
        } else {
          valid = pseudo_chol_solve(C, factor, centered_target, solution);
          if (valid) {
            const double alpha0 = arma::accu(solution) / t;
            W = solution.col(0) - v.col(0) * alpha0;
            alpha = target_mean + alpha0;
          }
        }
      }
    } else if (valid) {
      valid = pseudo_chol_solve(C, factor, target, solution);
      if (valid) W = solution.col(0);
    }
    valid = valid && std::isfinite(denominator) && denominator > 0 &&
      std::isfinite(alpha) && W.is_finite();
    long double score = 0;
    if (valid && !(intercept && constant)) {
      // Cancel lambda^2 algebraically before evaluating GCV.
      const long double ratio = static_cast<long double>(arma::norm(W, 2)) / denominator;
      score = n * ratio * ratio;
    }
    const double gcv = static_cast<double>(score);
    valid = valid && std::isfinite(gcv) && gcv >= 0;
    if (!valid) { skipped.push_back(lambda(j)); continue; }
    out.valid(j) = 1;
    out.scores(j) = gcv;
    // Strict comparison preserves supplied order for exact GCV ties.
    if (gcv < best) {
      best = gcv;
      out.W = W;
      out.alpha = alpha;
      out.lambda = lambda(j);
    }
  }
  if (!std::isfinite(out.lambda))
    Rcpp::stop("No supplied pseudo-direct ridge candidate has a valid positive-definite solve and GCV denominator");
  if (!skipped.empty()) {
    std::ostringstream message;
    message << "Skipped unusable pseudo-direct ridge candidate(s): ";
    for (std::size_t j = 0; j < skipped.size(); ++j) {
      if (j) message << ", ";
      message << skipped[j];
    }
    Rcpp::warning(message.str());
  }
  return out;
}

struct PseudoControls {
  double rho, eta, gamma, tau, btol, ftol, gtol;
  int maxitr, verbose, trial_budget;
};

void pseudo_check_controls(const PseudoControls& c) {
  if (!std::isfinite(c.rho) || c.rho <= 0 || c.rho >= 1 ||
      !std::isfinite(c.eta) || c.eta <= 0 || c.eta >= 1 ||
      !std::isfinite(c.gamma) || c.gamma < 0 || c.gamma > 1 ||
      !std::isfinite(c.tau) || c.tau <= 0 ||
      !std::isfinite(c.btol) || c.btol <= 0 ||
      !std::isfinite(c.ftol) || c.ftol <= 0 ||
      !std::isfinite(c.gtol) || c.gtol <= 0 || c.maxitr < 0 ||
      c.trial_budget < 0 || c.trial_budget > 30)
    Rcpp::stop("Invalid pseudo-direct optimizer controls: require 0 < rho, eta < 1, 0 <= gamma <= 1, positive finite steps and tolerances, and valid iteration budgets");
}

Rcpp::List pseudo_fit(arma::mat B, const arma::mat& X, const arma::vec& R,
                      const arma::vec& A, const arma::vec& grid,
                      const arma::vec& lambda, const double bw,
                      const bool intercept, const PseudoControls& control,
                      const bool collect, const int fail_gradient,
                      const int fail_trial) {
  pseudo_check_controls(control);
  const PseudoData data = pseudo_data(X, R, A, bw);
  PseudoEvaluation state = pseudo_evaluate(B, data, true);
  if (!state.valid) Rcpp::stop(state.message);
  if (pseudo_feasibility(B) > pseudo_feasibility_tolerance)
    Rcpp::stop("Initial pseudo-direct basis must be orthonormal within 1e-8");
  Rcpp::List trace;
  double Q = 1, reference = state.objective;
  double step = control.tau, previous_step = control.tau;
  int completed = 0;
  arma::mat changes(5, 2, arma::fill::zeros);
  std::string reason = "max_iterations";
  if (collect) trace.push_back(Rcpp::List::create(
    _["event"] = "initialize", _["B"] = B,
    _["objective"] = state.objective, _["gradient"] = state.gradient,
    _["gradient.norm"] = state.norm, _["reference"] = reference,
    _["reference.weight"] = Q, _["step"] = step,
    _["previous.step"] = previous_step));
  if (state.norm <= control.gtol) reason = "gradient_tolerance";
  while (completed < control.maxitr && reason == "max_iterations") {
    Rcpp::checkUserInterrupt();
    if (!std::isfinite(pseudo_feasibility(B)) ||
        pseudo_feasibility(B) > pseudo_feasibility_tolerance) {
      reason = "numerical_failure";
      break;
    }
    const PseudoCurve curve(B, state.gradient,
      !(B.n_cols < B.n_rows / 2));
    if (!std::isfinite(curve.slope) || curve.slope <= 0) {
      reason = "numerical_failure";
      break;
    }
    double trial_step = step;
    bool accepted = false;
    arma::mat accepted_B;
    PseudoEvaluation candidate;
    for (int trial = 0; trial < control.trial_budget; ++trial) {
      Rcpp::checkUserInterrupt();
      arma::mat trial_B;
      bool feasible = std::isfinite(trial_step) && trial_step > 0 &&
        curve.trial(B, trial_step, trial_B);
      PseudoEvaluation proposed;
      if (feasible) proposed = pseudo_evaluate(trial_B, data, false);
      if (completed + 1 == fail_trial) {
        proposed.valid = false;
        proposed.message = "Audit-injected trial failure";
      }
      const double threshold = reference - control.rho * trial_step * curve.slope;
      bool sufficient = feasible && proposed.valid && std::isfinite(threshold) &&
        proposed.objective <= threshold;
      bool derivative_valid = false;
      if (sufficient) {
        candidate = pseudo_evaluate(trial_B, data, true);
        if (completed + 1 == fail_gradient) {
          candidate.valid = false;
          candidate.message = "Audit-injected derivative failure";
        }
        derivative_valid = candidate.valid;
      }
      accepted = sufficient && derivative_valid;
      if (collect) trace.push_back(Rcpp::List::create(
        _["event"] = "trial", _["iteration"] = completed + 1,
        _["trial"] = trial + 1, _["B.start"] = B, _["B"] = trial_B,
        _["objective.start"] = state.objective, _["objective"] = proposed.objective,
        _["reference"] = reference, _["reference.weight"] = Q,
        _["step"] = trial_step, _["slope"] = curve.slope,
        _["threshold"] = threshold, _["feasible"] = feasible,
        _["valid"] = proposed.valid, _["derivative.valid"] = derivative_valid,
        _["accepted"] = accepted));
      if (accepted) { accepted_B = trial_B; break; }
      trial_step *= control.eta;
    }
    if (!accepted) {
      reason = "line_search_failed";
      if (collect) trace.push_back(Rcpp::List::create(
        _["event"] = "rollback", _["iteration"] = completed + 1,
        _["B"] = B, _["objective"] = state.objective,
        _["gradient"] = state.gradient, _["reference"] = reference,
        _["reference.weight"] = Q, _["step"] = step,
        _["previous.step"] = previous_step));
      break;
    }
    const arma::mat S = accepted_B - B;
    const arma::mat Y = candidate.residual - state.residual;
    const double B_change = arma::norm(S, "fro") / std::sqrt(static_cast<double>(B.n_rows));
    const double F_change = std::abs(candidate.objective - state.objective) /
      (std::abs(state.objective) + 1.0);
    bool valid_spectral;
    const double next_step = orthoDr::bb_step(arma::accu(S % S),
      std::abs(arma::accu(S % Y)), arma::accu(Y % Y),
      (completed + 1) % 2 == 0, previous_step, 1e10, valid_spectral);
    B = accepted_B;
    state = candidate;
    ++completed;
    step = previous_step = next_step;
    const double previous_Q = Q;
    Q = control.gamma * previous_Q + 1;
    reference = (control.gamma * previous_Q * reference + state.objective) / Q;
    changes((completed - 1) % 5, 0) = B_change;
    changes((completed - 1) % 5, 1) = F_change;
    if (collect) trace.push_back(Rcpp::List::create(
      _["event"] = "commit", _["iteration"] = completed, _["B"] = B,
      _["objective"] = state.objective, _["gradient"] = state.gradient,
      _["gradient.norm"] = state.norm, _["reference"] = reference,
      _["reference.weight"] = Q, _["step"] = step,
      _["previous.step"] = previous_step, _["accepted.step"] = trial_step,
      _["spectral.valid"] = valid_spectral));
    if (state.norm <= control.gtol) reason = "gradient_tolerance";
    // A single heavily backtracked step can be tiny while the next spectral
    // step is useful. Require small changes across five accepted updates.
    else if (completed >= 5 &&
      arma::mean(changes.col(0)) < control.btol &&
      arma::mean(changes.col(1)) < control.ftol) reason = "stagnation";
    if (control.verbose > 1 && completed % 10 == 0)
      Rcpp::Rcout << "At iteration " << completed << ", normalized MSE = " << state.objective << std::endl;
  }
  if (!std::isfinite(pseudo_feasibility(B)) ||
      pseudo_feasibility(B) > pseudo_feasibility_tolerance)
    reason = "numerical_failure";
  // Extraction is mandatory even after a numerical stop. Undefined targets or
  // ridge coefficients are errors, never a partially usable fit.
  const PseudoTargets targets = pseudo_targets(state, data, grid,
    pseudo_final_dose_factor, collect);
  const PseudoRidge ridge = pseudo_ridge(state.kernel, targets.dose, lambda, intercept);
  Rcpp::List out = Rcpp::List::create(
    _["B"] = B, _["W"] = ridge.W,
    _["fn"] = data.reward_scale * std::sqrt(state.objective),
    _["itr"] = completed, _["converge"] = reason == "gradient_tolerance",
    _["stop.reason"] = reason, _["gradient.norm"] = state.norm,
    _["objective.scaled"] = state.objective, _["reward.scale"] = data.reward_scale,
    _["lambda.init"] = ridge.lambda,
    _["objective.version"] = intercept ? "pseudo_mse_intercept_v1" : "pseudo_mse_v1",
    _["dose.bandwidth.scale"] = data.dose_scale * pseudo_final_dose_factor);
  if (intercept) out["intercept"] = ridge.alpha;
  if (collect) {
    out["trace"] = trace;
    out["audit.history"] = Rcpp::List::create(
      _["reference"] = reference, _["reference.weight"] = Q,
      _["step"] = step, _["previous.step"] = previous_step);
    out["grid.targets"] = targets.dose;
    out["grid.scores"] = targets.scores;
    out["GCV"] = ridge.scores;
    out["ridge.valid"] = ridge.valid;
  }
  if (control.verbose > 0)
    Rcpp::Rcout << "Pseudo-direct stop: " << reason << "; accepted updates: " << completed <<
      "; normalized MSE: " << state.objective << "; normalized gradient norm: " << state.norm << std::endl;
  return out;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List pdose_semi_solver(arma::mat B, const arma::mat& X,
                  const arma::colvec& R, const arma::colvec& A,
                  const arma::colvec a_seq, const arma::colvec lambda,
                  const double bw, double rho, double eta, double gamma,
                  double tau, double epsilon, double btol, double ftol,
                  double gtol, int maxitr, int verbose, int ncore,
                  bool intercept = false) {
  (void)epsilon; // Accepted for compatibility; gradients are analytic.
  (void)ncore;
  const PseudoControls control = {rho, eta, gamma, tau, btol, ftol, gtol,
    maxitr, verbose, 30};
  return pseudo_fit(B, X, R, A, a_seq, lambda, bw, intercept, control, false, -1, -1);
}

// Private fixed-state oracle access. Nonorthonormal B is allowed so coordinate
// central differences test the ambient derivative, not a reorthogonalization.
// [[Rcpp::export]]
Rcpp::List pdose_semi_evaluate(const arma::mat B, const arma::mat X,
                  const arma::colvec R, const arma::colvec A, double bw,
                  const arma::colvec a_seq, const arma::colvec lambda,
                  bool intercept = false, double cayley_step = 1e-4,
                  double dose_factor = 0.70710678118654752440) {
  const PseudoData data = pseudo_data(X, R, A, bw);
  const PseudoEvaluation state = pseudo_evaluate(B, data, true);
  if (!state.valid) Rcpp::stop(state.message);
  if (!std::isfinite(cayley_step) || cayley_step < 0)
    Rcpp::stop("Audit Cayley step must be finite and nonnegative");
  const PseudoTargets targets = pseudo_targets(state, data, a_seq, dose_factor, true);
  const PseudoRidge ridge = pseudo_ridge(state.kernel, targets.dose, lambda, intercept);
  const PseudoCurve dense(B, state.gradient, true), lowrank(B, state.gradient, false);
  arma::mat dense_B, lowrank_B;
  const bool dense_valid = dense.trial(B, cayley_step, dense_B);
  const bool lowrank_valid = lowrank.trial(B, cayley_step, lowrank_B);
  Rcpp::List out = Rcpp::List::create(
    _["objective.scaled"] = state.objective,
    _["fn"] = data.reward_scale * std::sqrt(state.objective),
    _["reward.scale"] = data.reward_scale, _["gradient"] = state.gradient,
    _["residual"] = state.residual, _["gradient.norm"] = state.norm,
    _["kernel"] = state.kernel, _["log.kernel"] = state.log_kernel,
    _["fitted.scaled"] = state.fitted, _["reward.scaled"] = data.q,
    _["grid.targets"] = targets.dose, _["grid.scores"] = targets.scores,
    _["W"] = ridge.W, _["lambda.init"] = ridge.lambda,
    _["GCV"] = ridge.scores, _["ridge.valid"] = ridge.valid,
    _["intercept"] = ridge.alpha,
    _["dose.bandwidth.scale"] = data.dose_scale * dose_factor);
  out["cayley.dense"] = Rcpp::List::create(_["slope"] = dense.slope,
    _["B"] = dense_B, _["valid"] = dense_valid);
  out["cayley.lowrank"] = Rcpp::List::create(_["slope"] = lowrank.slope,
    _["B"] = lowrank_B, _["valid"] = lowrank_valid);
  return out;
}

// Production and audit share the identical loop. Failure hooks only reject
// candidates and cannot install an invalid state.
// [[Rcpp::export]]
Rcpp::List pdose_semi_audit_solver(arma::mat B, const arma::mat& X,
                  const arma::colvec& R, const arma::colvec& A,
                  const arma::colvec a_seq, const arma::colvec lambda,
                  const double bw, double rho, double eta, double gamma,
                  double tau, double epsilon, double btol, double ftol,
                  double gtol, int maxitr, int verbose, int ncore,
                  bool intercept = false, bool trace = true,
                  int trial_budget = 30, int audit_fail_gradient_iteration = -1,
                  int audit_fail_trial_iteration = -1) {
  (void)epsilon;
  (void)ncore;
  const PseudoControls control = {rho, eta, gamma, tau, btol, ftol, gtol,
    maxitr, verbose, trial_budget};
  return pseudo_fit(B, X, R, A, a_seq, lambda, bw, intercept, control, trace,
    audit_fail_gradient_iteration, audit_fail_trial_iteration);
}
