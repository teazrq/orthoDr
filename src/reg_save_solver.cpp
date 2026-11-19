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
#include <limits>
#include "utilities.h"
#include "orthoDr_reg.h"
#include "reg_projection.h"

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

// Compute local deviations before taking moments in the supplied units.
static void save_covariance(const arma::mat& X, const arma::mat& kernel,
                            arma::mat& residual, arma::cube& covariance,
                            int ncore)
{
  const int N = X.n_rows;
  const arma::rowvec totals = sum(kernel, 0);

#pragma omp parallel for schedule(static) num_threads(ncore) if(ncore > 1)
  for (int i = 0; i < N; ++i) {
    arma::rowvec local_mean(X.n_cols, arma::fill::zeros);
    for (int j = 0; j < N; ++j) {
      if (kernel(i,j) > 0)
        local_mean += (X.row(j) - X.row(i)) * kernel(i,j);
    }
    local_mean /= totals(i);

    arma::rowvec deviation(X.n_cols);
    for (int j = 0; j < N; ++j) {
      if (kernel(i,j) > 0) {
        deviation = (X.row(j) - X.row(i)) - local_mean;
        covariance.slice(i) += deviation.t() * (deviation * kernel(i,j));
      }
    }
    covariance.slice(i) /= totals(i);
    residual.row(i) = -local_mean;
  }
}

static void save_response_moments(const arma::mat& X, const arma::mat& Y,
                                  arma::mat& residual, arma::cube& covariance,
                                  int ncore, double identity_weight,
                                  double covariance_weight)
{
  save_covariance(X, KernelDist_multi(Y, ncore, 1), residual, covariance, ncore);
  const arma::mat identity = identity_weight * arma::eye(X.n_cols, X.n_cols);
  covariance.each_slice([&](arma::mat& slice) { slice = identity - covariance_weight * slice; });
}

static arma::mat save_projection_kernel(const arma::mat& X, const arma::mat& B,
                                        double bw, int ncore)
{
  const arma::mat centered = X.each_row() - X.row(0);
  arma::mat projected = centered * B;
  const arma::rowvec scale = arma::stddev(projected, 0, 0) * bw * std::sqrt(2.0);
  const double eps = std::numeric_limits<double>::epsilon();
  const arma::rowvec roundoff = (eps * arma::max(arma::abs(centered), 0)) * arma::abs(B);
  if (!arma::any(roundoff > std::sqrt(eps) * scale)) {
    projected.each_row() /= scale;
    return ncore > 1 ? KernelDist_multi(projected, ncore, 1)
                     : KernelDist_single(projected, 1);
  }

  // A distant common origin can lose nearby differences at a narrow
  // bandwidth. Form each feature difference before projection in this case.
  return reg_pair_kernel(X, B, bw, ncore);
}

double save_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              const arma::mat& Exy,
              const arma::cube& Covxy,
              double bw,
              int ncore)
{
  // This function computes the estimation equations and its 2-norm for the semi-parametric dimensional reduction model

  int N = X.n_rows;
  int P = X.n_cols;
  arma::mat Ex(N, P, arma::fill::none);
  arma::cube Covxx(P, P, N, arma::fill::zeros);
  save_covariance(X, save_projection_kernel(X, B, bw, ncore), Ex, Covxx, ncore);

  arma::mat Est(P, P, arma::fill::zeros);

  for (int i=0; i < N; i++){
    Est += Covxy.slice(i) * (Exy.row(i).t()*Ex.row(i) - Covxx.slice(i));
  }

  return accu(pow(Est/N, 2));

}

// initial value

//  @title save_init
//  @name save_init
//  @description save initial value function
//  @keywords internal
// [[Rcpp::export]]
double save_init(const arma::mat& B,
                 const arma::mat& X,
                 const arma::mat& Y,
                 double bw,
                 int ncore,
                 double identity_weight = 1.0,
                 double covariance_weight = 1.0,
                 double moment_scale = 1.0)
{
  int N = X.n_rows;
  int P = B.n_rows;

  // initialize parallel computing

  checkCores(ncore, 0.0);

  arma::mat Exy(N, P, arma::fill::none);
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  save_response_moments(X, Y, Exy, Covxy, ncore, identity_weight, covariance_weight);

  // Initial function value

  double F = save_f(B, X, Y, Exy, Covxy, bw, ncore);

  return F * moment_scale;
}


//  @title semi-save solver \code{C++} function
//  @name save_solver
//  @description Solving the semi-save estimating equations. This is an internal function and should not be called directly.
//  @keywords internal
//  @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//  @param X A matrix of the parameters \code{X}
//  @param Y A matrix of the parameters \code{Y}
//  @param bw Kernel bandwidth for X
//  @param rho (don't change) Parameter for control the linear approximation in line search
//  @param eta (don't change) Factor for decreasing the step size in the backtracking line search
//  @param gamma (don't change) Parameter for updating C by Zhang and Hager (2004)
//  @param tau (don't change) Step size for updating
//  @param epsilon (don't change) Parameter for apprximating numerical gradient, if \code{g} is not given.
//  @param btol (don't change) The \code{$B$} parameter tolerance level
//  @param ftol (don't change) Functional value tolerance level
//  @param gtol (don't change) Gradient tolerance level
//  @param maxitr Maximum number of iterations
//  @param verbose Should information be displayed
//  @references Ma, Y. & Zhu, L. (2012). A semiparametric approach to dimension reduction. Journal of the American Statistical Association, 107(497), 168-179.
//  DOI: \url{https://dx.doi.org/10.1214\%2F12-AOS1072SUPP}.
//  @references Wen, Z. & Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434.
//  DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
// 
// [[Rcpp::export]]

List save_solver(arma::mat B,
                 arma::mat& X,
                 arma::mat& Y,
                 double bw,
                 double rho,
                 double eta,
                 double gamma,
                 double tau,
                 double epsilon,
                 double btol,
                 double ftol,
                 double gtol,
                 int maxitr,
                 int verbose,
                 int ncore,
                 double identity_weight = 1.0,
                 double covariance_weight = 1.0,
                 double moment_scale = 1.0)
{
  int N = X.n_rows;
  int P = B.n_rows;
  int ndr = B.n_cols;

  bool invH = true;

  if(ndr < P/2){
    invH = false;
  }

  // initialize parallel computing

  checkCores(ncore, verbose);

  arma::mat Exy(N, P, arma::fill::none);
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  save_response_moments(X, Y, Exy, Covxy, ncore, identity_weight, covariance_weight);

  // Build lambdas that capture the precalculated data and delegate to save_f

  auto f_eval = [&X, &Y, &Exy, &Covxy, bw, ncore, moment_scale](const arma::mat& B_) -> double {
    return save_f(B_, X, Y, Exy, Covxy, bw, ncore) * moment_scale;
  };

  auto g_eval = [&X, &Y, &Exy, &Covxy, bw, epsilon, ncore, moment_scale](arma::mat& B_,
                                                                       const double F0,
                                                                       arma::mat& G) -> void {
    numerical_gradient(B_, F0, G,
      [&X, &Y, &Exy, &Covxy, bw, moment_scale](const arma::mat& B_) -> double {
        return save_f(B_, X, Y, Exy, Covxy, bw, 1) * moment_scale;
      },
      epsilon, ncore);
  };

  // Initial function value and gradient

  double F = f_eval(B);

  arma::mat G(P, ndr, arma::fill::zeros);
  g_eval(B, F, G);

  // Delegate to the shared solver

  return orthoDr_solve(B, F, G, f_eval, g_eval,
                       P, ndr, maxitr,
                       rho, eta, gamma, tau,
                       btol, ftol, gtol,
                       verbose, invH);
}
