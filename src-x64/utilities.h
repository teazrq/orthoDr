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

#ifdef _OPENMP
#include <omp.h>
#define OMPMSG(...)
#else
#define omp_get_thread_num() 0
#define omp_get_max_threads() 1
#define OMPMSG(...) Rprintf("Package is not compiled with OpenMP (omp.h).\n")
#endif

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

#ifndef orthoDr_utility
#define orthoDr_utility

double dmax(double a, double b);
double imax(int a, int b);
double dmin(double a, double b);
double imin(int a, int b);

void checkCores(int& ncore, int verbose);

arma::mat KernelDist_multi(const arma::mat& X, int ncore, double diag);
arma::mat KernelDist_single(const arma::mat& X, double diag);
arma::mat EpanKernelDist_multi(const arma::mat& X, int ncore, double diag);
arma::mat EpanKernelDist_single(const arma::mat& X, double diag);

NumericMatrix KernelDist_cross(const arma::mat& TestX, const arma::mat& X);

#endif
