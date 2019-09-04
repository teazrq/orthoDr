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

#ifndef orthoDr_reg
#define orthoDr_reg

double local_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              double bw,
              int ncore);
			  
void local_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Y,
            double bw,
            double epsilon,
            int ncore);

Rcpp::List local_solver(arma::mat B,
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
                int ncore);

double phd_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              const arma::cube& XX,
              double bw,
              int ncore);

void phd_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Y,
            const arma::cube& XX,
            double bw,
            double epsilon,
            int ncore);

double phd_init(const arma::mat& B,
                 const arma::mat& X,
                 const arma::mat& Y,
                 double bw,
                 int ncore);

Rcpp::List phd_solver(arma::mat B,
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
                 int ncore);

double save_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              const arma::mat& Exy,
              const arma::cube& Covxy,
              double bw,
              int ncore);

void save_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Y,
            const arma::mat& Exy,
            const arma::cube& Covxy,
            double bw,
            double epsilon,
            int ncore);

double save_init(const arma::mat& B,
                 const arma::mat& X,
                 const arma::mat& Y,
                 double bw,
                 int ncore);
				 
Rcpp::List save_solver(arma::mat B,
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
                 int ncore);
				 
double seff_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              const arma::mat& kernel_matrix_y,
              double bw,
              int ncore);

void seff_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Y,
            const arma::mat& kernel_matrix_y,
            double bw,
            double epsilon,
            int ncore);

double seff_init(const arma::mat& B,
                const arma::mat& X,
                const arma::mat& Y,
                double bw,
                int ncore);

Rcpp::List seff_solver(arma::mat B,
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
                int ncore);

double sir_f(const arma::mat& B,
             const arma::mat& X,
             const arma::mat& Exy,
             double bw,
             int ncore);

void sir_g(arma::mat& B,
           const double F0,
           arma::mat& G,
           const arma::mat& X,
           const arma::mat& Exy,
           double bw,
           double epsilon,
           int ncore);

double sir_init(const arma::mat& B,
                const arma::mat& X,
                const arma::mat& Y,
                double bw,
                int ncore);

Rcpp::List sir_solver(arma::mat B,
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
                int ncore);

				
#endif
