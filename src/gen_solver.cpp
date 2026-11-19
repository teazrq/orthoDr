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
#include "orthoDr_gen.h"
#include "solver_numerics.h"

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

double gen_f(const arma::mat& B, Rcpp::Function f, Environment env)
{
  Rcpp::RObject x = Rcpp_eval(f(B), env);
  if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || Rf_xlength(x) != 1)
    Rcpp::stop("fn must return a single number at every evaluation");
  double value = Rcpp::as<double>(x);
  if (!std::isfinite(value)) Rcpp::stop("fn must return a finite number at every evaluation");
  return value;
}


void gen_g(const arma::mat& B, arma::mat& G, Rcpp::Function g, Environment env)
{
  Rcpp::RObject x = Rcpp_eval(g(B), env);

  int P = B.n_rows;
  int ndr = B.n_cols;

  if (TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP)
    Rcpp::stop("grad must return a numeric matrix at every evaluation");
  SEXP dims = Rf_getAttrib(x, R_DimSymbol);
  if (Rf_length(dims) != 2 || INTEGER(dims)[0] != P || INTEGER(dims)[1] != ndr)
    Rcpp::stop("grad must return a matrix with the same dimension as B at every evaluation");
  G = Rcpp::as<arma::mat>(x);
  if (!G.is_finite()) Rcpp::stop("grad must return a finite matrix at every evaluation");

  return;
}


void gen_g_approx(const arma::mat& B, arma::mat& G, Rcpp::Function f, Environment env, double epsilon)
{
  int P = B.n_rows;
  int ndr = B.n_cols;

  arma::mat NewB = B;

  for (int j = 0; j < ndr; j++)
  {
    for(int i = 0; i < P; i++)
    {
      NewB(i,j) = B(i,j) + epsilon;
      const double above = gen_f(NewB, f, env);
      NewB(i,j) = B(i,j) - epsilon;
      G(i,j) = (above - gen_f(NewB, f, env)) / (2 * epsilon);

      // reset from original
      NewB(i,j) = B(i,j);
    }
  }

  return;
}

arma::mat tangent_gradient(const arma::mat& B, const arma::mat& G)
{
  const arma::mat BtG = B.t() * G;
  arma::mat projected = G - B * (0.5 * BtG + 0.5 * BtG.t());
  // Remove any residual normal component left by subtraction of a large term.
  const arma::mat residual = B.t() * projected;
  return projected - B * (0.5 * residual + 0.5 * residual.t());
}

//  @title General solver \code{C++} function
//  @name gen_solver
//  @description A general purpose optimization solver with orthogonality constraint. For details, please see the original \code{MATLAB} code by Wen and Yin (2013). This is an internal function and should not be called directly.
//  @keywords internal
//  @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//  @param f A function that calculates the objective function value. The first argument should be \code{B}. Returns a single value.
//  @param g A function that calculates the gradient. The first argument should be \code{B}. Returns a matrix with the same dimension as \code{B}. If not specified, then the numerical approximation is used.
//  @param env Environment passed to the Rcpp function for evaluating \code{f} and \code{g}
//  @param useg If true, the gradient is calculated using \code{g} function, otherwise numerically approximated
//  @param rho (don't change) Parameter for control the linear approximation in line search
//  @param eta (don't change) Factor for decreasing the step size in the backtracking line search
//  @param gamma (don't change) Parameter for updating C by Zhang and Hager (2004)
//  @param tau (don't change) Step size for updating
//  @param epsilon (don't change) Parameter for approximating numerical gradient, if \code{g} is not given.
//  @param btol (don't change) The \code{$B$} parameter tolerance level
//  @param ftol (don't change) Functional value tolerance level
//  @param gtol (don't change) Gradient tolerance level
//  @param maxitr Maximum number of iterations
//  @param verbose Should information be displayed
//  @return A list with components: \code{B} (optimal matrix), \code{fn} (final objective value), \code{fn_Seq} (sequence of objective values), \code{itr} (number of iterations), \code{converge} (convergence indicator).
//  @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//  @references Zhang, H. and Hager, W. W., "A nonmonotone line search technique and its application to unconstrained optimization." SIAM J. Optim. 14 (2004): 1043–1056. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//  @examples
//  # This function should be called internally. When having all objects pre-computed, one can call
//  # gen_solver(B, f, g, env, useg, rho, eta, gamma, tau, epsilon, btol, ftol, gtol, maxitr, verbose)
//  # to solve for the parameters B.
// [[Rcpp::export]]


List gen_solver(arma::mat B,
                 Rcpp::Function f,
                 Rcpp::Function g,
                 Environment env,
                 int useg,
                 double rho,
                 double eta,
                 double gamma,
                 double tau,
                 double epsilon,
                 double btol,
                 double ftol,
                 double gtol,
                 int maxitr,
                 int verbose)
{

  const int P = B.n_rows;
  const int ndr = B.n_cols;
  const bool dense = ndr >= P / 2;
  if (B.is_empty() || !B.is_finite() ||
      orthoDr::feasibility_error(B) > orthoDr::feasibility_tolerance)
    Rcpp::stop("Initial B must have orthonormal columns");

  auto objective = [&](const arma::mat& b) { return gen_f(b, f, env); };
  auto projected_gradient = [&](const arma::mat& b, double value) {
    arma::mat G(P, ndr);
    if (useg)
      gen_g(b, G, g, env);
    else
      gen_g_approx(b, G, f, env, epsilon);
    orthoDr::check_state(value, G);
    return tangent_gradient(b, G);
  };
  auto gradient_norm = [](const arma::mat& D) {
    const double norm = arma::norm(D, "fro");
    if (!std::isfinite(norm)) Rcpp::stop("Constrained gradient norm must be finite");
    return norm;
  };

  double F = objective(B);
  arma::vec F_seq(maxitr, arma::fill::zeros);
  arma::mat changes(5, 2, arma::fill::zeros);
  arma::mat D = projected_gradient(B, F);
  double nrmG = gradient_norm(D);
  double best_norm = nrmG;
  int completed = 0;
  int last_improvement = 0;
  const char* reason = "iteration limit reached";
  double Q = 1;
  double Cval = F;

  if (verbose > 1)
    Rcout << "Initial value,   F = " << F << std::endl;

  while (nrmG > gtol && completed < maxitr) {
    const arma::mat BP = B;
    const double FP = F;
    // Removing the normal component leaves the skew generator unchanged.
    // Scaling keeps the reduced system well scaled even for large gradients.
    const arma::mat unit = D / nrmG;
    const arma::mat direction = unit - B * (unit.t() * B);
    const double slope = rho * nrmG * arma::accu(unit % direction);
    const orthoDr::CayleyCurve curve(B, unit, dense);
    double scaled_step = tau * nrmG;
    if (!orthoDr::backtrack(B, F, scaled_step, eta, Cval, slope, curve, objective)) {
      reason = "line search failed to find a feasible sufficient-decrease step";
      break;
    }

    const arma::mat nextD = projected_gradient(B, F);
    const arma::mat S = B - BP;
    const arma::mat Y = nextD - D;
    const double accepted_tau = scaled_step / nrmG;
    D = nextD;
    nrmG = gradient_norm(D);
    F_seq[completed] = F;
    changes(completed % 5, 0) = arma::norm(S, "fro") / std::sqrt(static_cast<double>(P));
    changes(completed % 5, 1) = std::abs(FP - F) / (std::abs(FP) + 1);
    ++completed;
    if (verbose > 1 && completed % 10 == 0)
      Rcout << "At iteration " << completed << ", F = " << F << std::endl;
    if (nrmG <= gtol) break;
    if (nrmG < best_norm) {
      best_norm = nrmG;
      last_improvement = completed;
    }
    if (completed - last_improvement >= 5 && arma::mean(changes.col(0)) < btol &&
        arma::mean(changes.col(1)) < ftol) {
      reason = "parameter and objective changes stagnated";
      break;
    }

    bool valid_bb;
    tau = orthoDr::bb_step(arma::accu(S % S), std::abs(arma::accu(S % Y)),
                           arma::accu(Y % Y), completed % 2 == 0,
                           accepted_tau > 0 ? accepted_tau : tau, 1e20, valid_bb);
    Q = gamma * Q + 1;
    Cval += (F - Cval) / Q;
  }

  const bool converged = nrmG <= gtol;
  if (verbose > 0){
    if (!converged) Rcout << reason << std::endl;
    Rcout << "number of iterations: " << completed << std::endl;
    Rcout << "norm of functional value: " << F << std::endl;
    Rcout << "norm of gradient: " << nrmG << std::endl;
    Rcout << "norm of feasibility: " << orthoDr::feasibility_error(B) << std::endl;
  }

  List ret;
  ret["B"] = B;
  ret["fn"] = F;
  ret["fn_Seq"] = F_seq;
  ret["itr"] = completed;
  ret["converge"] = converged;
  return (ret);
}
