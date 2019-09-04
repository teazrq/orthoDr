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

#include "utilities.h"

#ifndef orthoDr_pdose
#define orthoDr_pdose

double pdose_direct_f(const arma::mat& B,
                   const arma::mat& X,
                   const arma::colvec& A,
                   const arma::colvec& a_seq,
                   const arma::colvec& R,
                   const double bw,
                   const arma::colvec& W,
                   int ncore);
				   
				   
void pdose_direct_w(const arma::mat& B,
                 const arma::mat& X,
                 const arma::colvec& A,
                 const arma::colvec& a_seq,
                 const arma::mat& a_dist,
                 const arma::colvec& R,
                 const double bw,
                 arma::colvec& W,
                 const arma::colvec& lambda,
                 int ncore);
				 
void pdose_direct_g(arma::mat& B,
                 const double F0,
                 arma::mat& G,
                 const arma::mat& X,
                 const arma::colvec& A,
                 const arma::colvec& a_seq,
                 const arma::colvec& R,
                 const double bw,
                 const arma::colvec& W,
                 const double lambda0,
                 const double epsilon,
                 int ncore);

Rcpp::List pdose_direct_solver(arma::mat B,
                      const arma::mat X,
                      const arma::colvec A,
                      const arma::mat a_dist,
                      const arma::colvec a_seq,
                      const arma::colvec R,
                      const arma::colvec lambda,
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
                      int ncore);

arma::vec dosepred(arma::mat B,
                   arma::mat X,
                   arma::mat X_test,
                   double bw,
                   arma::colvec W);
				   
double pdose_semi_f(const arma::mat& B,
                      const arma::mat& X,
                      const arma::colvec& R,
                      const arma::colvec& A,
                      const double bw,
                      int ncore);

void  pdose_semi_g( const arma::mat& B,
                 const double F0,
                 arma::mat& G,
                 const arma::mat& X,
                 const arma::colvec& R,
                 const arma::colvec& A,
                 const double bw,
                 const double epsilon,
                 int ncore);

Rcpp::List pdose_semi_solver(arma::mat& B,
                  const arma::mat& X,
                  const arma::colvec& R,
                  const arma::colvec& A,
                  const arma::mat a_dist,
                  const arma::colvec a_seq,
                  const arma::colvec lambda,
                  const double bw,
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
                  int ncore);				 
				   
#endif
