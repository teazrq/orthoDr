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
};

struct DirectEvaluation {
  bool valid;
  std::string message;
  double objective;
  double value;
  double penalty;
  double norm_B;
  double norm_W;
  arma::vec dose;
  arma::mat gradient_B;
  arma::vec gradient_W;
  arma::vec curvature_dose;
  DirectKernel geometry;

  DirectEvaluation() : valid(false), objective(NA_REAL), value(NA_REAL),
    penalty(NA_REAL), norm_B(NA_REAL), norm_W(NA_REAL) {}
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
  return out.kernel.is_finite();
}

DirectEvaluation direct_evaluate(const arma::mat& B, const arma::mat& X,
                                 const arma::vec& A, const arma::vec& R,
                                 const double bw, const arma::vec& W,
                                 const double lambda_reward,
                                 const bool gradients = true) {
  DirectEvaluation out;
  const arma::uword n = X.n_rows;
  if (A.n_elem != n || R.n_elem != n || W.n_elem != n ||
      !A.is_finite() || !R.is_finite() || !W.is_finite() ||
      !std::isfinite(lambda_reward) || lambda_reward < 0) {
    out.message = "Direct-dose outcomes and coefficients must be finite with matching lengths, and lambda.reward must be nonnegative";
    return out;
  }
  if (!direct_kernel(B, X, bw, out.geometry, out.message)) return out;
  const double dose_scale = arma::stddev(A) * bw;
  if (!std::isfinite(dose_scale) || dose_scale <= 0) {
    out.message = "Observed doses must have a positive finite bandwidth scale";
    return out;
  }
  out.dose = out.geometry.kernel * W;
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
    out.gradient_W = out.geometry.kernel.t() * gradient_dose +
      2.0 * (lambda_reward * W);
    // Each column i has dose d_i = sum_j K(j,i) W_j. Both kernel
    // endpoints and the sample standard deviation move with B.
    const arma::mat sensitivity = log_sensitivity +
      out.geometry.kernel % (W * gradient_dose.t());
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
    if (!out.gradient_B.is_finite() || !out.gradient_W.is_finite() ||
        !std::isfinite(out.norm_B) || !std::isfinite(out.norm_W)) {
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
  arma::mat weighted_kernel = state.geometry.kernel;
  weighted_kernel.each_col() %= state.curvature_dose;
  hessian = state.geometry.kernel.t() * weighted_kernel;
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
                                double lambda_reward, double cayley_step = 0) {
  const DirectEvaluation state = direct_evaluate(B, X, A, R, bw, W, lambda_reward);
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
  }
  return out;
}

namespace {

struct DirectInitialization {
  arma::vec W;
  arma::vec targets;
  arma::vec gcv;
  double lambda;
  int count;
};

DirectInitialization direct_initialize(const DirectKernel& geometry,
                                       const arma::vec& A,
                                       const arma::vec& grid,
                                       const arma::vec& R,
                                       const arma::vec& lambda,
                                       const double bw) {
  DirectInitialization out;
  const arma::uword n = A.n_elem;
  out.count = 1;
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
                  const arma::vec& direction_W = arma::vec()) {
  if (!collect) return;
  trace.push_back(Rcpp::List::create(
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
    Rcpp::Named("direction.W") = direction_W));
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
                      const Rcpp::Nullable<Rcpp::NumericVector>& W_initial) {
  direct_validate(B, X, A, R, a_dist, grid, lambda, bw, lambda_reward, control);
  DirectKernel initial_geometry;
  std::string message;
  if (!direct_kernel(B, X, bw, initial_geometry, message)) Rcpp::stop(message);
  DirectInitialization initialization;
  if (W_initial.isNotNull()) {
    initialization.W = Rcpp::as<arma::vec>(W_initial.get());
    initialization.lambda = lambda(0);
    initialization.count = 0;
  } else {
    initialization = direct_initialize(initial_geometry, A, grid, R, lambda, bw);
  }
  arma::vec W = initialization.W;
  DirectEvaluation state = direct_evaluate(B, X, A, R, bw, W, lambda_reward);
  if (!state.valid) Rcpp::stop("Initial direct-dose evaluation failed: " + state.message);
  Rcpp::List trace;
  if (collect) {
    direct_trace(trace, true, "initialize", "", 0, B, W, B, W,
                  state.objective, state, NA_REAL, 0, true);
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
    if (state.norm_B <= control.gtol && state.norm_W <= control.gtol) {
      reason = "gradient_tolerance";
      break;
    }
    const arma::mat previous_B = B;
    const arma::vec previous_W = W;
    const double previous_objective = state.objective;
    arma::mat proposed_B = B;
    arma::vec proposed_W = W;
    DirectEvaluation proposed = state;
    bool failed = false;
    std::string failed_block;

    if (proposed.norm_W > control.gtol) {
      failed_block = "W";
      const arma::vec gradient = proposed.gradient_W;
      std::string direction_type;
      const arma::vec direction = direct_direction_W(proposed, lambda_reward, direction_type);
      const double slope = arma::dot(gradient, direction);
      if (!std::isfinite(slope) || slope >= 0) {
        reason = std::isfinite(slope) ? "stagnation" : "numerical_failure";
        failed = true;
      } else {
        bool accepted = false;
        double step = W_step;
        for (int trial = 0; trial < 30; ++trial) {
          Rcpp::checkUserInterrupt();
          const arma::vec trial_W = W + step * direction;
          DirectEvaluation candidate;
          if (trial_W.is_finite())
            candidate = direct_evaluate(B, X, A, R, bw, trial_W, lambda_reward, false);
          else candidate.message = "Coefficient trial is nonfinite";
          const bool moved = trial_W.n_elem == W.n_elem && arma::any(trial_W != W);
          const double threshold = proposed.objective + control.rho * step * slope;
          accepted = moved && candidate.valid && std::isfinite(threshold) &&
            candidate.objective <= threshold;
          direct_trace(trace, collect, "trial", "W", completed + 1, B, W,
                        B, trial_W, proposed.objective, candidate, slope, step, accepted,
                        "", direction_type, direction);
          if (candidate.valid && !moved) {
            reason = "stagnation";
            failed = true;
            break;
          }
          if (accepted) {
            proposed_W = trial_W;
            proposed = candidate;
            W_step = direct_next_step(step, control.eta, 1.0);
            break;
          }
          step *= control.eta;
        }
        if (!accepted && !failed) {
          reason = "w_line_search_failed";
          failed = true;
        }
      }
    } else {
      direct_trace(trace, collect, "skip", "W", completed + 1, B, W, B, W,
                    proposed.objective, proposed, 0, 0, true);
    }

    // Recompute the direction derivative at provisional coefficients. A
    // derivative failure rolls back their accepted objective-only trial.
    if (!failed) {
      proposed = direct_evaluate(proposed_B, X, A, R, bw, proposed_W, lambda_reward);
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
            candidate = direct_evaluate(trial_B, X, A, R, bw, proposed_W, lambda_reward, false);
          else candidate.message = "Cayley trial is nonfinite or unsolved";
          const bool moved = feasible && arma::accu(trial_B != B) > 0;
          const double threshold = proposed.objective + control.rho * step * slope;
          accepted = moved && candidate.valid && std::isfinite(threshold) &&
            candidate.objective <= threshold;
          direct_trace(trace, collect, "trial", "B", completed + 1, B, proposed_W,
                        trial_B, proposed_W, proposed.objective, candidate, slope, step, accepted);
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
                    B, proposed_W, proposed.objective, proposed, 0, 0, true);
    }

    // Both fresh gradients and all components belong to the final pair,
    // including at the last budgeted sweep, before the atomic commit.
    if (!failed) {
      proposed = direct_evaluate(proposed_B, X, A, R, bw, proposed_W, lambda_reward);
      if (!proposed.valid) {
        reason = "numerical_failure";
        failed_block = "B";
        failed = true;
      }
    }
    if (failed) {
      direct_trace(trace, collect, "rollback", failed_block, completed + 1,
                    proposed_B, proposed_W, B, W, previous_objective, state,
                    NA_REAL, 0, false, reason);
      break;
    }
    B = proposed_B;
    W = proposed_W;
    state = proposed;
    ++completed;
    direct_trace(trace, collect, "commit", "", completed, previous_B,
                  previous_W, B, W, previous_objective, state, NA_REAL, 0, true);
    if (state.norm_B <= control.gtol && state.norm_W <= control.gtol) {
      reason = "gradient_tolerance";
      break;
    }
    const double B_change = arma::norm(B - previous_B, "fro") /
      std::max(1.0, arma::norm(previous_B, "fro"));
    const double W_change = arma::norm(W - previous_W, 2) /
      std::max(1.0, arma::norm(previous_W, 2));
    const double objective_change = std::abs(state.objective - previous_objective) /
      (1.0 + std::abs(previous_objective));
    if (B_change <= control.btol && W_change <= control.btol &&
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
  if (state.norm_B <= control.gtol && state.norm_W <= control.gtol)
    reason = "gradient_tolerance";
  const DirectEvaluation fresh = direct_evaluate(B, X, A, R, bw, W, lambda_reward);
  if (fresh.valid) state = fresh;
  else {
    reason = "numerical_failure";
    state.norm_B = NA_REAL;
    state.norm_W = NA_REAL;
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
    Rcpp::Named("objective.version") = "direct_joint_v1");
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
                              int verbose, int ncore, double lambda_reward = 0.1) {
  // gamma and epsilon are retained only for internal calling compatibility.
  (void)gamma;
  (void)epsilon;
  (void)ncore;
  const DirectControls control = {rho, eta, tau, btol, ftol, gtol,
                                  maxitr, verbose, 30};
  return direct_fit(B, X, A, a_dist, a_seq, R, lambda, bw, lambda_reward,
                     control, false, R_NilValue);
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
                                    Rcpp::Nullable<Rcpp::NumericVector> W_initial = R_NilValue) {
  (void)gamma;
  (void)epsilon;
  (void)ncore;
  const DirectControls control = {rho, eta, tau, btol, ftol, gtol,
                                  maxitr, verbose, B_trial_budget};
  return direct_fit(B, X, A, a_dist, a_seq, R, lambda, bw, lambda_reward,
                     control, trace, W_initial);
}
