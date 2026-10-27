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

#include "utilities.h"

#ifndef orthoDr_pdose
#define orthoDr_pdose

Rcpp::List pdose_direct_evaluate(const arma::mat B,
                                const arma::mat X,
                                const arma::colvec A,
                                const arma::colvec R,
                                double bw,
                                const arma::colvec W,
                                double lambda_reward,
                                double cayley_step,
                                bool intercept,
                                Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds,
                                double alpha);

Rcpp::List pdose_direct_policy_direction(const arma::vec gradient,
                                        const arma::mat hessian,
                                        double alpha,
                                        const Rcpp::NumericVector intercept_bounds);

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
                      int ncore,
                      double lambda_reward,
                      bool intercept,
                      Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds);

Rcpp::List pdose_direct_audit_solver(arma::mat B,
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
                                    int ncore,
                                    double lambda_reward,
                                    bool trace,
                                    int B_trial_budget,
                                    Rcpp::Nullable<Rcpp::NumericVector> W_initial,
                                    bool intercept,
                                    Rcpp::Nullable<Rcpp::NumericVector> intercept_bounds,
                                    Rcpp::Nullable<Rcpp::NumericVector> alpha_initial);

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
