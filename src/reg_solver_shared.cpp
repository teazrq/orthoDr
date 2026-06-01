//    ----------------------------------------------------------------
//
//    Orthogonality Constrained Optimization for Dimension Reduction
//    (orthoDr)
//
//    Shared solver infrastructure for regression methods.
//    Contains the Wen-Yin curvilinear search loop and numerical
//    gradient computation used by all 5 regression solvers.
//
//    ----------------------------------------------------------------

#include "orthoDr_reg.h"
#include <functional>
#include <algorithm>

// Numerical gradient via finite differences (thread-safe)
// Each thread gets its own copy of B; writes to unique G(i,j) elements
// via collapse(2) schedule(static) — no data races.
void numerical_gradient(arma::mat& B,
                        const double F0,
                        arma::mat& G,
                        const std::function<double(const arma::mat&)>& f_eval,
                        double epsilon,
                        int ncore)
{
  int P = B.n_rows;
  int ndr = B.n_cols;

#pragma omp parallel num_threads(ncore)
{
  // thread-local copy of B
  arma::mat NewB(P, ndr);
  NewB = B;

  #pragma omp for collapse(2) schedule(static)
  for (int j = 0; j < ndr; j++)
  for (int i = 0; i < P; i++)
  {
    double temp = B(i, j);
    NewB(i, j) = B(i, j) + epsilon;

    // finite difference gradient
    G(i,j) = (f_eval(NewB) - F0) / epsilon;

    // reset
    NewB(i, j) = temp;
  }
}
}

// Wen-Yin curvilinear search solver for orthogonality-constrained problems
// (Wen & Yin 2013, Algorithm 1)
//
// Thread safety: The main loop is sequential. f_eval and g_eval may use
// OpenMP internally but must not write to shared state outside their
// own parallel regions.
Rcpp::List orthoDr_solve(arma::mat& B,
                         double F,
                         arma::mat& G,
                         const std::function<double(const arma::mat&)>& f_eval,
                         const std::function<void(arma::mat&, const double, arma::mat&)>& g_eval,
                         int P,
                         int ndr,
                         int maxitr,
                         double rho,
                         double eta,
                         double gamma,
                         double tau,
                         double btol,
                         double ftol,
                         double gtol,
                         int verbose,
                         bool invH)
{
  // Guard against invalid maxitr (R validates, but defend at C++ level too)
  maxitr = std::max(5, maxitr);

  arma::mat crit(maxitr, 3, arma::fill::zeros);

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
  int itr;
  arma::mat BP;
  double FP;
  arma::mat GP;
  arma::mat dtXP;
  arma::mat diag_n(P, P);
  arma::mat aa;
  arma::mat S;
  double BDiff;
  double FDiff;
  arma::mat Y_Y;
  double SY;

  arma::mat eye2P(2*ndr, 2*ndr);
  eye2P.eye();

  if (verbose > 1)
    Rcpp::Rcout << "Initial value,   F = " << F << std::endl;

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

      F = f_eval(B);
      g_eval(B, F, G);

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

    dtX = G - B * GX;
    nrmG = norm(dtX, "fro");

    S = B - BP;
    BDiff = norm(S, "fro")/sqrt((double) P);
    FDiff = std::abs(FP - F)/(std::abs(FP)+1);

    Y_Y = dtX - dtXP;
    SY = std::abs(accu(S % Y_Y));

    if(itr%2 == 0){
      tau = accu(S % S)/SY;
    }else{
      tau = SY/accu(Y_Y % Y_Y);
    }

    tau = std::max(std::min(tau, 1e10), 1e-20);
    crit(itr-1,0) = nrmG;
    crit(itr-1,1) = BDiff;
    crit(itr-1,2) = FDiff;

    if (verbose > 1 && (itr % 10 == 0) )
      Rcpp::Rcout << "At iteration " << itr << ", F = " << F << std::endl;

    if (itr >= 5)
    {
      arma::mat mcrit(5, 3);
      for (int i=0; i<5; i++)
      {
        mcrit.row(i) = crit.row(itr-i-1);
      }

      if ( (BDiff < btol && FDiff < ftol) || (nrmG < gtol) || ((mean(mcrit.col(1)) < btol) && (mean(mcrit.col(2)) < ftol)) )
      {
        if (verbose > 0) Rcpp::Rcout << "converge" << std::endl;
        break;
      }
    }

    double Qp = Q;
    Q = gamma * Qp + 1;
    Cval = (gamma*Qp*Cval + F)/Q;

  }

  if(itr>=maxitr){
    if (verbose > 0) Rcpp::Rcout << "exceed max iteration before convergence ... " << std::endl;
  }

  arma::mat diag_P(ndr,ndr);
  diag_P.eye();
  double feasi = norm(B.t() * B - diag_P, "fro");

  if (verbose > 0){
    Rcpp::Rcout << "number of iterations: " << itr << std::endl;
    Rcpp::Rcout << "norm of functional value: " << F << std::endl;
    Rcpp::Rcout << "norm of gradient: " << nrmG << std::endl;
    Rcpp::Rcout << "norm of feasibility: " << feasi << std::endl;
  }

  Rcpp::List ret;
  ret["B"] = B;
  ret["fn"] = F;
  ret["itr"] = itr;
  ret["converge"] = (itr<maxitr);
  return (ret);
}
