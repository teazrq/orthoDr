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

double semi_pt_f(const arma::mat& B,
                      const arma::mat& X,
                      const arma::colvec& R,
                      const arma::colvec& A,
                      const double bw,
                      int ncore)

{
  int N = X.n_rows;
  int ndr = B.n_cols;

  arma::mat BX  = X * B;
  BX.insert_cols(1, A);

  arma::rowvec BX_scale = stddev(BX, 0, 0)*bw*sqrt(2);

  for (int j=0; j<ndr+1; j++)
    BX.col(j) /= BX_scale(j);

  arma::mat kernel_matrix;

  if (ncore > 1)
    kernel_matrix=  KernelDist_multi(BX,ncore, 1);
  else
    kernel_matrix =  KernelDist_single(BX,1);

  arma::vec kernel_mean(N);
  for (int i = 0; i < N; i++) {
    kernel_mean(i) = sum(R % kernel_matrix.col(i))/sum(kernel_matrix.col(i));
  }

  double mse = sqrt(mean( (R-kernel_mean) % (R-kernel_mean)) );

  return mse;
}



void  semi_pt_g( const arma::mat& B,
                 const double F0,
                 arma::mat& G,
                 const arma::mat& X,
                 const arma::colvec& R,
                 const arma::colvec& A,
                 const double bw,
                 const double epsilon,
                 int ncore)
{
  // This function computes the gradiant of the estimation equations

  int P = B.n_rows;
  int ndr = B.n_cols;

#pragma omp parallel num_threads(ncore)
{
  // create one copy of B for each thread
  arma::mat NewB(P, ndr);
  NewB = B;

#pragma omp for collapse(2) schedule(static)
  for (int j = 0; j < ndr; j++)
    for(int i = 0; i < P; i++)
    {
      // small increment
      double temp = B(i, j);
      NewB(i, j) = B(i, j) + epsilon;

      // calculate gradiant
      G(i,j) = (semi_pt_f(NewB, X, R, A, bw, ncore) - F0) / epsilon;

      // reset
      NewB(i,j) = temp;
    }
}

return;
}




//' @title pseudo_pt_solver
//' @name pseudo_pt_solver
//' @description The main optimization function for personalized dose finding in the dimensional reduction framewrok, the Pseudo-direct Learning method.
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param X The covariate matrix
//' @param R The perosnalzied medicine reward
//' @param A observed dose levels
//' @param bw A Kernel bandwidth, assuming each variable have unit variance
//' @param rho (don't change) Parameter for control the linear approximation in line search
//' @param eta (don't change) Factor for decreasing the step size in the backtracking line search
//' @param gamma (don't change) Parameter for updating C by Zhang and Hager (2004)
//' @param tau (don't change) Step size for updating
//' @param epsilon (don't change) Parameter for approximating numerical gradient
//' @param btol (don't change) The \code{$B$} parameter tolerance level
//' @param ftol (don't change) Estimation equation 2-norm tolerance level
//' @param gtol (don't change) Gradient tolerance level
//' @param maxitr Maximum number of iterations
//' @param verbose Should information be displayed
//' @return The optimizer \code{B} for the esitmating equation.
//' @references Zhou, W., Zhu, R. "A Parsimonious Personalized Dose Model vis Dimension Reduction." (2018)  \url{https://arxiv.org/abs/1802.06156}.
//' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
// [[Rcpp::export]]

List semi_pt_solver(arma::mat& B,
                  const arma::mat& X,
                  const arma::colvec& R,
                  const arma::colvec& A,
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
                  int ncore)
{

  int P = B.n_rows;
  int ndr = B.n_cols;

  arma::mat crit(maxitr,3);
  bool invH = true;
  arma::mat eye2P(2*ndr,2*ndr);

  if(ndr < P/2){
    invH = false;
    eye2P.eye();
  }

  if (ncore > 1) OMPMSG(1);
  int haveCore = omp_get_max_threads();
  if (ncore <= 0) ncore = haveCore;

  if (ncore > haveCore)
  {
    if (verbose) Rcout << "Do not have " <<  ncore << " cores, use maximum " << haveCore << " cores." << std::endl;
    ncore = haveCore;
  }

  // Initial function value and gradient, prepare for iterations

  double F = semi_pt_f(B, X, R, A, bw, ncore);

  arma::mat G(P, ndr);
  G.fill(0);
  semi_pt_g(B, F, G, X, R, A, bw, epsilon, ncore);

  //return G;

  arma::mat GX = G.t() * B;
  arma::mat GXT;
  arma::mat H;
  arma::mat RX;
  arma::mat U;
  arma::mat V;
  arma::mat VU;
  arma::mat VX;

  if(invH){
    GXT = G * B.t();
    H = 0.5 * (GXT - GXT.t());
    RX = H * B;
  }else{
    U = join_rows(G, B);
    V = join_rows(B, -G);
    VU = V.t() * U;
    VX = V.t() * B;
  }

  arma::mat dtX = G - B * GX;
  double nrmG = norm(dtX, "fro");

  double Q = 1;
  double Cval = F;

  // main iteration
  int itr = 1;
  arma::mat BP;
  double FP;
  arma::mat GP;
  arma::mat dtXP;
  arma::mat diag_n(P, P);
  arma::mat aa;
  arma::mat S;
  double BDiff;
  double FDiff;
  arma::mat Y;
  double SY;

  if (verbose > 1)

  for(itr = 1; itr < maxitr + 1; itr++){
    BP = B;
    FP = F;
    GP = G;
    dtXP = dtX;

    int nls = 1;
    double deriv = rho * nrmG * nrmG;

    while(true){
      if(invH){
        diag_n.eye();
        B = solve(diag_n + tau * H, BP - tau * RX);
      }else{
        aa = solve(eye2P + 0.5 * tau * VU, VX);
        B = BP - U * (tau * aa);
      }

      F = semi_pt_f(B, X, R, A, bw, ncore);
      semi_pt_g(B, F, G, X, R, A, bw, epsilon, ncore);

      if((F <= (Cval - tau*deriv)) || (nls >= 5)){
        break;
      }
      tau = eta * tau;
      nls = nls + 1;
    }

    GX = G.t() * B;

    if(invH){
      GXT = G * B.t();
      H = 0.5 * (GXT - GXT.t());
      RX = H * B;
    }else{
      U = join_rows(G, B);
      V = join_rows(B, -G);
      VU = V.t() * U;
      VX = V.t() * B;
    }

    dtX = G - B * GX; // GX, dtX, nrmG slightly different from those of R code
    nrmG = norm(dtX, "fro");

    S = B - BP;
    BDiff = norm(S, "fro")/sqrt((double) P);
    FDiff = std::abs(FP - F)/(std::abs(FP)+1);

    Y = dtX - dtXP;
    SY = std::abs(accu(S % Y));

    if(itr%2 == 0){
      tau = accu(S % S)/SY;
    }else{
      tau = SY/accu(Y % Y);
    }

    tau = dmax(dmin(tau, 1e10), 1e-20);
    crit(itr-1,0) = nrmG;
    crit(itr-1,1) = BDiff;
    crit(itr-1,2) = FDiff;

    if (verbose > 1 && (itr % 10 == 0) )
      Rcout << "At iteration " << itr << ", F = " << F << std::endl;

    if (itr >= 5) // so I will run at least 5 iterations before checking for convergence
    {
      arma::mat mcrit(5, 3);
      for (int i=0; i<5; i++)
      {
        mcrit.row(i) = crit.row(itr-i-1);
      }

      if ( (BDiff < btol && FDiff < ftol) || (nrmG < gtol) || ((mean(mcrit.col(1)) < btol) && (mean(mcrit.col(2)) < ftol)) )
      {
        if (verbose > 0) Rcout << "converge" << std::endl;
        break;
      }
    }

    double Qp = Q;
    Q = gamma * Qp + 1;
    Cval = (gamma*Qp*Cval + F)/Q;

  }

	if(itr>=maxitr){
		Rcout << "exceed max iteration before convergence ... " << std::endl;
  }

  arma::mat diag_P(ndr,ndr);
  diag_P.eye();
  double feasi = norm(B.t() * B - diag_P, "fro");

  if (verbose > 0){
    Rcout << "number of iterations: " << itr << std::endl;
    Rcout << "functional value: " << std::setprecision(6) << F << std::endl;
    Rcout << "norm of gradient: " << nrmG << std::endl;
    Rcout << "norm of feasibility: " << feasi << std::endl;
  }

  List ret;
  ret["B"] = B;
  ret["fn"] = F;
  ret["itr"] = itr;
  ret["converge"] = (itr<maxitr);
  return (ret);
}
