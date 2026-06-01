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

  // E[X | BX]

  arma::mat Ex(N,P, arma::fill::zeros);

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Ex.row(i) += X.row(j)*kernel_matrix_x(i,j);
    }
    Ex.row(i) /= Kx(i);
  }

  // cov[X | XB]

  arma::cube Covxx(P, P, N, arma::fill::zeros);

  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Covxx.slice(i) += (X.row(j).t() * X.row(j))*kernel_matrix_x(i,j);
    }
    Covxx.slice(i) /= Kx(i);
    Covxx.slice(i) -= Ex.row(i).t() * Ex.row(i);
    Ex.row(i) = X.row(i) - Ex.row(i);
  }

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
                 int ncore)
{
  int N = X.n_rows;
  int P = B.n_rows;

  // initialize parallel computing

  checkCores(ncore, 0.0);

  //precalculate

  arma::mat kernel_matrix_y = KernelDist_multi(Y, ncore, 1);

  arma::rowvec Ky = sum(kernel_matrix_y, 0);

  // X - E[X | Y]
  arma::mat Exy(N, P, arma::fill::zeros);

  // I - cov[X | Y]
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  arma::mat diag = arma::eye(P, P);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Exy.row(i) += X.row(j)*kernel_matrix_y(i,j);
      Covxy.slice(i) += (X.row(j).t() * X.row(j)) * kernel_matrix_y(i,j);
    }

    // E[X | Y]
    Exy.row(i) /= Ky(i);

    // E[XX | Y]
    Covxy.slice(i) /= Ky(i);

    // I - cov[X | Y]
    Covxy.slice(i) = diag - Covxy.slice(i) + Exy.row(i).t()*Exy.row(i);

    // X - E[X | Y]
    Exy.row(i) = X.row(i) - Exy.row(i);
  }

  // Initial function value

  double F = save_f(B, X, Y, Exy, Covxy, bw, ncore);

  return F;
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
                 int ncore)
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

  //precalculate

  arma::mat kernel_matrix_y = KernelDist_multi(Y, ncore, 1);

  arma::rowvec Ky = sum(kernel_matrix_y, 0);

  // X - E[X | Y]
  arma::mat Exy(N, P, arma::fill::zeros);

  // I - cov[X | Y]
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  arma::mat diag = arma::eye(P, P);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Exy.row(i) += X.row(j)*kernel_matrix_y(i,j);
      Covxy.slice(i) += (X.row(j).t() * X.row(j)) * kernel_matrix_y(i,j);
    }

    // E[X | Y]
    Exy.row(i) /= Ky(i);

    // E[XX | Y]
    Covxy.slice(i) /= Ky(i);

    // I - cov[X | Y]
    Covxy.slice(i) = diag - Covxy.slice(i) + Exy.row(i).t()*Exy.row(i);

    // X - E[X | Y]
    Exy.row(i) = X.row(i) - Exy.row(i);
  }

  // Build lambdas that capture the precalculated data and delegate to save_f

  auto f_eval = [&X, &Y, &Exy, &Covxy, bw, ncore](const arma::mat& B_) -> double {
    return save_f(B_, X, Y, Exy, Covxy, bw, ncore);
  };

  auto g_eval = [&X, &Y, &Exy, &Covxy, bw, epsilon, ncore](arma::mat& B_,
                                                             const double F0,
                                                             arma::mat& G) -> void {
    numerical_gradient(B_, F0, G,
      [&X, &Y, &Exy, &Covxy, bw](const arma::mat& B_) -> double {
        return save_f(B_, X, Y, Exy, Covxy, bw, 1);
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
