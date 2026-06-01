//    ----------------------------------------------------------------
//
//    Semiparametric efficient estimator (Ma & Zhu 2013, Section 3.2).
//
//    Two-step procedure:
//      Step 1: Obtain root-n consistent initial β̃ (SAVE plug-in)
//      Steps 2-3: Estimate η̂₂, ∂η̂₂/∂u, and Ê using data at β̃ᵀx (ONCE)
//      Step 4: Optimize S_eff — re-evaluate η̂₂ and Ê at each candidate
//              βᵀx using kernel interpolation with FIXED β̃ᵀx_j.
//
//    Key: the nonparametric functions are ESTIMATED from β̃ᵀx_j,
//    but EVALUATED at the current βᵀx_i during optimization.
//
//    ----------------------------------------------------------------

#include <RcppArmadillo.h>
#include "utilities.h"
#include "orthoDr_reg.h"

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

// Efficient score objective (Step 4).
// B_tilde_X_ref, Y_ref, Ky_ref, bx_scale_ref = pre-computed at β̃.
// At each candidate B, re-evaluates η̂₂ and Ê at βᵀx_i via interpolation.
double seff_f_fixed(const arma::mat& B,
                     const arma::mat& X,
                     const arma::mat& B_tilde_X_ref,   // β̃ᵀX_j (fixed)
                     const arma::mat& Y_ref,            // Y_j (fixed)
                     const arma::mat& Ky_ref,            // Y-kernel (fixed)
                     double bw,
                     double bw_y,                        // Y-kernel bandwidth
                     int ncore)
{
  int N = X.n_rows;
  int P = X.n_cols;
  int ndr = B.n_cols;

  // Current projection
  arma::mat BX = X * B;
  arma::rowvec BX_scale = stddev(BX, 0, 0) * bw * sqrt(2.0);
  for (int j = 0; j < ndr; j++) BX.col(j) /= BX_scale(j);

  // Scale B_tilde_X_ref the same way as BX for fair comparison
  arma::mat B_ref_scaled = B_tilde_X_ref;
  arma::rowvec B_ref_scale = stddev(B_ref_scaled, 0, 0) * bw * sqrt(2.0);
  for (int j = 0; j < ndr; j++) B_ref_scaled.col(j) /= B_ref_scale(j);

  // Kernel interpolation weights: K(BX_i, B_tilde_X_j) for all i,j
  // This is an N×N matrix
  arma::mat Kx_interp(N, N, arma::fill::zeros);
  
  // Efficient kernel computation (only need lower triangle + diagonal)
  #pragma omp parallel for schedule(static) num_threads(ncore)
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double d2 = 0;
      for (int k = 0; k < ndr; k++) {
        double d = BX(i, k) - B_ref_scaled(j, k);
        d2 += d * d;
      }
      Kx_interp(i, j) = exp(-d2);
    }
  }

  // Ê[X|βᵀX_i] via Nadaraya-Watson using β̃ᵀx_j as reference
  arma::mat Ex(N, P, arma::fill::zeros);
  #pragma omp parallel for schedule(static) num_threads(ncore)
  for (int i = 0; i < N; i++) {
    double w_sum = 0;
    for (int j = 0; j < N; j++) {
      double w = Kx_interp(i, j);
      w_sum += w;
      for (int l = 0; l < P; l++)
        Ex(i, l) += X(j, l) * w;
    }
    for (int l = 0; l < P; l++)
      Ex(i, l) /= w_sum;
  }

  // ∂log η̂₂(Y_i, βᵀX_i) via local linear regression
  // Regress Ky[:,i] on (B_tilde_X_j - BX_i) with weights K(BX_i, B_tilde_X_j)
  arma::mat score_hat(N, ndr);
  arma::vec sqK_col(N);
  arma::mat X_w(N, ndr + 1);

  for (int i = 0; i < N; i++) {
    sqK_col = sqrt(Kx_interp.row(i).t());  // sqrt weights for this i
    
    // Design matrix columns
    X_w.col(0) = sqK_col;
    for (int k = 0; k < ndr; k++) {
      for (int j = 0; j < N; j++)
        X_w(j, k + 1) = (B_ref_scaled(j, k) - BX(i, k)) * sqK_col(j);
    }

    // Weighted response
    arma::mat resp = Ky_ref.col(i) % sqK_col;

    arma::mat beta_hat = arma::solve(X_w, resp, arma::solve_opts::fast);
    double a_i = beta_hat(0, 0);
    for (int k = 0; k < ndr; k++)
      score_hat(i, k) = beta_hat(k + 1, 0) / a_i;
  }

  // Efficient score: Σ (X_i − Ê[X_i|βᵀX_i]) × ∂log η̂₂(Y_i, βᵀX_i)/∂u
  arma::mat Est(P, ndr, arma::fill::zeros);
  for (int i = 0; i < N; i++) {
    Est += (X.row(i) - Ex.row(i)).t() * score_hat.row(i);
  }

  return accu(pow(Est / N, 2));
}

// initial value
// [[Rcpp::export]]
double seff_init(const arma::mat& B,
                const arma::mat& X,
                const arma::mat& Y,
                double bw,
                int ncore)
{
  int N = X.n_rows;
  int ndr = B.n_cols;
  checkCores(ncore, 0.0);

  // Pre-compute reference quantities at initial B
  arma::mat B_ref = X * B;
  
  // Y-kernel bandwidth: Ma & Zhu (2013, Section 3.2) says b = N^{-1/7} * sd(Y)
  // Since Y is pre-scaled in R (sd = N^{1/(ndr+5)}/sqrt(2)), the effective
  // bandwidth in Y_R units is: bw_y = N^{1/(ndr+5) - 1/7}
  double bw_y = std::pow((double)N, 1.0 / (ndr + 5.0) - 1.0 / 7.0);
  arma::mat Y_scaled = Y / bw_y;
  arma::mat Ky(N, N);
  if (ncore > 1)
    Ky = KernelDist_multi(Y_scaled, ncore, 1);
  else
    Ky = KernelDist_single(Y_scaled, 1);

  double F = seff_f_fixed(B, X, B_ref, Y, Ky, bw, bw_y, ncore);
  return F;
}


// Efficient semiparametric solver (Ma & Zhu 2013, Section 3.2).
// [[Rcpp::export]]
List seff_solver(arma::mat B,           // initial β̃ (root-n consistent)
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
                int ncore)
{
  int N = X.n_rows;
  int P = B.n_rows;
  int ndr = B.n_cols;

  bool invH = true;
  if (ndr < P / 2) invH = false;

  checkCores(ncore, verbose);

  // Step 2-3: Pre-compute reference quantities at initial β̃
  arma::mat B_tilde_X = X * B;  // fixed: β̃ᵀX_j for all j
  
  // Y-kernel bandwidth: Ma & Zhu (2013, Section 3.2) says b = N^{-1/7} * sd(Y)
  // Since Y is pre-scaled in R (sd = N^{1/(ndr+5)}/sqrt(2)), the effective
  // bandwidth in Y_R units is: bw_y = N^{1/(ndr+5) - 1/7}
  double bw_y = std::pow((double)N, 1.0 / (ndr + 5.0) - 1.0 / 7.0);
  arma::mat Y_scaled = Y / bw_y;
  arma::mat Ky(N, N);
  if (ncore > 1)
    Ky = KernelDist_multi(Y_scaled, ncore, 1);
  else
    Ky = KernelDist_single(Y_scaled, 1);

  // Step 4: Optimize with fixed reference, re-evaluating at each candidate B
  auto f_eval = [&X, &B_tilde_X, &Y, &Ky, bw, bw_y, ncore](const arma::mat& B_) -> double {
    return seff_f_fixed(B_, X, B_tilde_X, Y, Ky, bw, bw_y, ncore);
  };

  auto g_eval = [&f_eval, epsilon, ncore](arma::mat& B_,
                                           const double F0,
                                           arma::mat& G) -> void {
    numerical_gradient(B_, F0, G, f_eval, epsilon, ncore);
  };

  double F = f_eval(B);
  arma::mat G(P, ndr, arma::fill::zeros);
  g_eval(B, F, G);

  return orthoDr_solve(B, F, G, f_eval, g_eval,
                       P, ndr, maxitr,
                       rho, eta, gamma, tau,
                       btol, ftol, gtol,
                       verbose, invH);
}
