#include <RcppArmadillo.h>

using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]

//' @title The prediction function for the personalized direct learning dose model
//' @name Dosepred
//' @description Predict the fitted dose from the direct learning dose model
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param X The covariate matrix
//' @param X_test The test covariate matrix
//' @param bw A Kernel bandwidth, assuming each variable have unit variance
//' @param w The kernel ridge regression coefficient 
//' @return The predicted dose
// [[Rcpp::export]]
arma::vec Dosepred(arma::mat B,
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

  arma::rowvec BX_scale = stddev(BX, 0, 0)*bw*sqrt(2);
  arma::rowvec BX_test_scale = stddev(BX, 0, 0)*bw*sqrt(2);

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



