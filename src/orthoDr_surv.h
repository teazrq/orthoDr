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

#ifndef orthoDr_surv
#define orthoDr_surv

double surv_dm_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Phit,
              const arma::vec& Fail_Ind,
              double bw,
              int ncore);
			  
void surv_dm_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Phit,
            const arma::vec& Fail_Ind,
            double bw,
            const double epsilon,
            int ncore);
			
Rcpp::List surv_dm_solver(arma::mat B,
                 const arma::mat& X,
                 const arma::mat& Phit,
                 const arma::vec& Fail_Ind,
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
				 
double surv_dn_f(const arma::mat& B,
                const arma::mat& X,
                const arma::mat& Phit,
                const arma::vec& Fail_Ind,
                double bw,
                int ncore);

void surv_dn_g(arma::mat& B,
            double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Phit,
            const arma::vec& Fail_Ind,
            double bw,
            double epsilon,
            int ncore);

Rcpp::List surv_dn_solver(arma::mat B,
                 const arma::mat& X,
                 const arma::mat& Phit,
                 const arma::vec& Fail_Ind,
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

double surv_forward_f(const arma::mat& B,
                const arma::mat& X,
                const arma::vec& Fail_Ind,
                double bw,
                int ncore);

void surv_forward_g(arma::mat& B,
            double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::vec& Fail_Ind,
            double bw,
            double epsilon,
            int ncore);

Rcpp::List surv_forward_solver(arma::mat B,
                 const arma::mat& X,
                 const arma::vec& Fail_Ind,
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

				 
#endif
