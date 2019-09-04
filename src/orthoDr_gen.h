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

#ifndef orthoDr_gen
#define orthoDr_gen

double gen_f(arma::mat &B, Rcpp::Function f, Rcpp::Environment env);
void gen_g(arma::mat B, arma::mat &G, Rcpp::Function g, Rcpp::Environment env);
void gen_g_approx(arma::mat &B, arma::mat &G, Rcpp::Function f, Rcpp::Function g, Rcpp::Environment env, double epsilon);

Rcpp::List gen_solver(arma::mat B,
                 Rcpp::Function f,
                 Rcpp::Function g,
                 Rcpp::Environment env,
                 int useg,
                 double rho,
                 double eta,
                 double gamma,
                 double tau,
                 double epsilon,
                 double btol,
                 double ftol,
                 double gtol,
                 int maxitr,
                 int verbose);
				 
#endif
