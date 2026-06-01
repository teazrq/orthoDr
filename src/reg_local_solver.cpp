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
#include "orthoDr_reg.h"

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

// local method f value function
//  @keywords internal
// [[Rcpp::export]]
double local_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              double bw,
              int ncore)
{
  int N = X.n_rows;
  int P = X.n_cols;
  int ndr = B.n_cols;

  arma::mat BX = X * B;

  arma::rowvec BX_scale = stddev(BX, 0, 0)*bw*sqrt(2.0);
  for (int j=0; j<ndr; j++)
    BX.col(j) /= BX_scale(j);

  arma::mat kernel_matrix_x(N, N);

  if (ncore > 1)
    kernel_matrix_x = KernelDist_multi(BX, ncore, 1);
  else
    kernel_matrix_x = KernelDist_single(BX, 1);

  arma::rowvec Kx = sum(kernel_matrix_x, 0);

  // E[X | BX] via Nadaraya-Watson kernel regression
  arma::mat Ex(N, P, arma::fill::zeros);

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Ex.row(i) += X.row(j)*kernel_matrix_x(i,j);
    }
    Ex.row(i) /= Kx(i);
  }

  // E[Y | BX] and dE[Y]/d(BX) via local linear regression of Y on BX
  // For each i, solve: min_{a_i, b_i} sum_j { Y_j - a_i - b_i'(BX_j - BX_i) }^2 K(BX_i, BX_j)
  // This gives a_i = E_n[Y | BX_i] and b_i = dE_n[Y | BX] / d(BX) |_{BX_i}
  //
  // Ref: Ma & Zhu (2013) Annals of Statistics 41(1), Section 3, Eq (5)-(6)

  arma::vec a(N);
  arma::mat b(N, ndr);
  arma::mat X_w(N, ndr + 1);

  // use sqrt(kernel) for weighted least squares
  arma::mat sqK = sqrt(kernel_matrix_x);

  for(int i=0; i<N; i++){

    for(int j=0; j<N; j++)
      for (int k=0; k<ndr; k++)
        X_w(j, k+1) = BX(j, k) - BX(i, k);

    X_w.col(0) = sqK.col(i);

    for (int k=1; k<ndr+1; k++)
      X_w.col(k) = X_w.col(k) % sqK.col(i);

    arma::mat beta_hat = arma::solve(X_w, Y % sqK.col(i), arma::solve_opts::fast);

    a(i) = beta_hat(0, 0);

    for(int k=0; k<ndr; k++)
      b(i, k) = beta_hat(k+1, 0);
  }

  // Estimating equation: Psi(B) = (1/N) sum_i {X_i - E[X|BX_i]} {Y_i - a_i} b_i'
  // Objective: f(B) = ||Psi(B)||_F^2 / N^2

  double ret = 0;
  arma::mat Est(P, ndr, arma::fill::zeros);

  for(int i=0; i<N; i++){
    Est += (X.row(i) - Ex.row(i)).t() * (Y(i, 0) - a(i)) * b.row(i);
  }

  ret = accu(pow(Est, 2)) / N / N;

  return ret;
}

// Solve local semiparametric estimating equations (Ma & Zhu 2013)
// Uses Wen-Yin curvilinear search for orthogonality constraints
// [[Rcpp::export]]

List local_solver(arma::mat B,
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
                int ncore)
{
  int P = B.n_rows;
  int ndr = B.n_cols;

  bool invH = true;
  if(ndr < P/2){
    invH = false;
  }

  // initialize parallel computing
  checkCores(ncore, verbose);

  // Bind parameters into lambdas for the shared solver infrastructure
  std::function<double(const arma::mat&)> f_eval =
    [&X, &Y, bw, ncore](const arma::mat& B_eval) -> double {
      return local_f(B_eval, X, Y, bw, ncore);
    };

  std::function<void(arma::mat&, const double, arma::mat&)> g_eval =
    [&f_eval, epsilon, ncore](arma::mat& B_eval, const double F0, arma::mat& G) {
      numerical_gradient(B_eval, F0, G, f_eval, epsilon, ncore);
    };

  // Initial function value and gradient
  double F = f_eval(B);
  arma::mat G(P, ndr, arma::fill::zeros);
  g_eval(B, F, G);

  // Delegate to shared Wen-Yin solver
  return orthoDr_solve(B, F, G, f_eval, g_eval, P, ndr, maxitr, rho, eta,
                       gamma, tau, btol, ftol, gtol, verbose, invH);
}
