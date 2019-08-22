//    ----------------------------------------------------------------
//
//    Orthogonality Constrained Optimization for Dimension Reduction
//    (orthoDr)
//
//    This program is free software; you can redistribute it and/or
//    modify it under the terms of the GNU General Public License
//    as published by the Free Software Foundation; either version 3
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

using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' @title The prediction function for the personalized direct learning dose model
//' @name dosepred
//' @description Predict the fitted dose from the direct learning dose model
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param X The covariate matrix
//' @param X_test The test covariate matrix
//' @param bw A Kernel bandwidth, assuming each variable have unit variance
//' @param w The kernel ridge regression coefficient
//' @return The predicted dose
// [[Rcpp::export]]

arma::vec dosepred(arma::mat B,
                   arma::mat X,
                   arma::mat X_test,
                   double bw,
                   arma::colvec W)

{
  int N = X.n_rows;
  int N_test = X_test.n_rows;
  int ndr = B.n_cols;

  arma::mat BX = X * B;
  arma::mat BX_test = X_test * B;

  arma::rowvec BX_scale = stddev(BX, 0, 0)*bw*sqrt(2.0);
  arma::rowvec BX_test_scale = stddev(BX, 0, 0)*bw*sqrt(2.0);

  for (int j=0; j<ndr; j++)
    BX.col(j) /= BX_scale(j);

  for (int j=0; j<ndr; j++)
    BX_test.col(j) /= BX_test_scale(j);

  // calculate distance

  arma::mat kernel_matrix_X(N, N_test);

  for (int i = 0; i < N_test; i++)
  {
    for (int j = 0; j < N; j ++)
    {
      kernel_matrix_X(j,i) = exp(-sum(pow(BX.row(j)-BX_test.row(i),2)));
    }
  }

  arma::colvec Dose;
  Dose = kernel_matrix_X.t() * W;

  return Dose;

}



