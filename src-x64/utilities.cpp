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

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double dmax(double a, double b)
{
  if (a > b)
    return a;

  return b;
}

double imax(int a, int b)
{
  if (a > b)
    return a;

  return b;
}

double dmin(double a, double b)
{
  if (a > b)
    return b;

  return a;
}

double imin(int a, int b)
{
  if (a > b)
    return b;

  return a;
}

// check cores

void checkCores(int& ncore, int verbose)
{
  if (ncore > 1) OMPMSG(1);
  int haveCore = omp_get_max_threads();
  if (ncore <= 0) ncore = haveCore;

  if (ncore > haveCore)
  {
    if (verbose) Rcout << "Do not have " <<  ncore << " cores, use maximum " << haveCore << " cores." << std::endl;
    ncore = haveCore;
  }
}

// kernel distance functions

arma::mat KernelDist_multi(const arma::mat& X, int ncore, double diag)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for (int i = 0; i < (int) ceil((double) N /2); i++)
  {
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      kernel_matrix(j,i) = exp(-sum(pow(X.row(i)-X.row(j),2)));
      kernel_matrix(i,j) = kernel_matrix(j,i);
    }

    int m = N - i - 1;
    kernel_matrix(m, m) = diag;
    for (int j = 0; j < m; j ++)
    {
      kernel_matrix(j,m) = exp(-sum(pow(X.row(m)-X.row(j),2)));
      kernel_matrix(m,j) = kernel_matrix(j,m);
    }
  }

  return(kernel_matrix);
}


arma::mat KernelDist_single(const arma::mat& X, double diag)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N);

  for (int i = 0; i < N; i++)
  {
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      kernel_matrix(j,i) = exp(-sum(pow(X.row(i)-X.row(j),2)));
      kernel_matrix(i,j) = kernel_matrix(j,i);
    }
  }

  return(kernel_matrix);
}


arma::mat EpanKernelDist_single(const arma::mat& X, double diag)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N, arma::fill::zeros);
  double u;

  for (int i = 0; i < N; i++)
  {
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      u = sum(pow(X.row(i)-X.row(j),2));

      if (u > -1 && u < 1)
        kernel_matrix(j,i) = pow((1-u*u), 3);

      kernel_matrix(i,j) = kernel_matrix(j,i);
    }
  }
  return(kernel_matrix);
}



arma::mat EpanKernelDist_multi(const arma::mat& X, int ncore, double diag)
{
  int N = X.n_rows;
  arma::mat kernel_matrix(N, N, arma::fill::zeros);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for (int i = 0; i < (int) ceil((double) N /2); i++)
  {
    double u;
    kernel_matrix(i, i) = diag;
    for (int j = 0; j < i; j ++)
    {
      u = sum(pow(X.row(i)-X.row(j),2));

      if (u > -1 && u < 1)
        kernel_matrix(j,i) = pow((1-u*u), 3);

      kernel_matrix(i,j) = kernel_matrix(j,i);
    }

    int m = N - i - 1;
    kernel_matrix(m, m) = diag;
    for (int j = 0; j < m; j ++)
    {
      u = sum(pow(X.row(i)-X.row(j),2));

      if (u > -1 && u < 1)
        kernel_matrix(j,i) = pow((1-u*u), 3);

      kernel_matrix(m,j) = kernel_matrix(j,m);
    }
  }

  return(kernel_matrix);
}


//' @title KernelDist_cross
//' @name KernelDist_cross
//' @description Calculate the kernel distance between testing data and training data
//' @keywords internal
//' @param TestX testing data
//' @param X training data
// [[Rcpp::export]]
NumericMatrix KernelDist_cross(const arma::mat& TestX, const arma::mat& X)
{

  int N = X.n_rows;
  int TestN = TestX.n_rows;

  NumericMatrix kernel_matrix(TestN, N);

  for (int i = 0; i < TestN; i++){
  for (int j = 0; j < N; j++)
  {
    kernel_matrix(i,j) = exp(-sum(pow(TestX.row(i)-X.row(j),2)));
  }
  }
  return(kernel_matrix);
}
