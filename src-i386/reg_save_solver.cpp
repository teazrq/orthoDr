#include <RcppArmadillo.h>
#include "utilities.h"

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

double save_f(const arma::mat& B,
              const arma::mat& X,
              const arma::mat& Y,
              const arma::mat& Exy,
              const arma::cube& Covxy,
              double bw,
              int ncore)
{
  // This function computes the estimation equations and its 2-norm for the semi-parametric dimensional reduction model
  // OptCpp1(B, G, X, Phit_cpp, inRisk, kernel.bw.scale, Fail.Ind)

  int N = X.n_rows;
  int P = X.n_cols;
  int ndr = B.n_cols;

  arma::mat BX = X * B;

  arma::rowvec BX_scale = stddev(BX, 0, 0)*bw*sqrt(2.0);

  for (int j=0; j<ndr; j++)
    BX.col(j) /= BX_scale(j);

 // Rcout << stddev(BX, 0, 0) << std::endl;

  //Rcout << "Start" << std::endl;

  arma::mat kernel_matrix_x(N, N);

  if (ncore > 1)
    kernel_matrix_x = KernelDist_multi(BX, ncore, 1);
  else
    kernel_matrix_x = KernelDist_single(BX, 1);

  arma::rowvec Kx = sum(kernel_matrix_x, 0);

  //Rcout << "After Kx" << std::endl;

  // E[X | BX]

  arma::mat Ex(N,P, arma::fill::zeros);

// #pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Ex.row(i) += X.row(j)*kernel_matrix_x(i,j);
    }
    Ex.row(i) /= Kx(i);
    // Ex.row(i) = X.row(i) - Ex.row(i);
  }

  // arma::rowvec EEx = sum(Ex, 0);

  // cov[X | XB]

  arma::cube Covxx(P, P, N, arma::fill::zeros);

// #pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Covxx.slice(i) += (X.row(j).t() * X.row(j))*kernel_matrix_x(i,j);
    }
    Covxx.slice(i) /= Kx(i);
    Covxx.slice(i) -= Ex.row(i).t() * Ex.row(i);
    Ex.row(i) = X.row(i) - Ex.row(i);
  }

  //Rcout << Ex << std::endl;

  arma::mat Est(P, P, arma::fill::zeros);

  for (int i=0; i < N; i++){
    Est += Covxy.slice(i) * (Exy.row(i).t()*Ex.row(i) - Covxx.slice(i));

    //Rcout << Covxy.slice(i) * (Exy.row(i).t()*Ex.row(i) - Covxx.slice(i)) << std::endl;
  }

/*
  arma::cube EE(P, P, N, arma::fill::zeros);

  // #pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      EE.slice(i) += (X.row(j).t() * X.row(j))*kernel_matrix_x(i,j);
    }
    EE.slice(i) /= Kx(i);

    EE.slice(i) = (X.row(i).t() * X.row(i) - EE.slice(i));
  }

  arma::mat Est(P, P, arma::fill::zeros);

  for (int i=0; i < N; i++){
    Est += EE.slice(i);
  }
*/

  return accu(pow(Est/N, 2));

}

void save_g(arma::mat& B,
            const double F0,
            arma::mat& G,
            const arma::mat& X,
            const arma::mat& Y,
            const arma::mat& Exy,
            const arma::cube& Covxy,
            double bw,
            double epsilon,
            int ncore)
  {
  int P = B.n_rows;
  int ndr = B.n_cols;

#pragma omp parallel num_threads(ncore)
{
  // create one copy of B for each thread
  arma::mat NewB(P, ndr);
  NewB = B;

  #pragma omp for collapse(2) schedule(static)
  for (int j = 0; j < ndr; j++)
    for (int i = 0; i < P; i++)
    {  // small increment

      double temp = B(i, j);
      NewB(i, j) = B(i, j) + epsilon;

      // calculate gradiant
      G(i,j) = (save_f(NewB, X, Y, Exy, Covxy, bw, 1) - F0) / epsilon;

      // reset
      NewB(i, j) = temp;
    }
}
  return;
}

// initial value

//' @title save_init
//' @name save_init
//' @description save initial value function
//' @keywords internal
// [[Rcpp::export]]
double save_init(const arma::mat& B,
                 const arma::mat& X,
                 const arma::mat& Y,
                 double bw,
                 int ncore)
{
  int N = X.n_rows;
  int P = B.n_rows;

  // initialize parallel computing

  checkCores(ncore, 0.0);

  //precalculate

  arma::mat kernel_matrix_y = KernelDist_multi(Y, ncore, 1);

  arma::rowvec Ky = sum(kernel_matrix_y, 0);

  // X - E[X | Y]
  arma::mat Exy(N, P, arma::fill::zeros);

  // I - cov[X | Y]
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  arma::mat diag = arma::eye(P, P);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Exy.row(i) += X.row(j)*kernel_matrix_y(i,j);
      Covxy.slice(i) += (X.row(j).t() * X.row(j)) * kernel_matrix_y(i,j);
    }

    // E[X | Y]
    Exy.row(i) /= Ky(i);

    // E[XX | Y]
    Covxy.slice(i) /= Ky(i);

    // I - cov[X | Y]
    Covxy.slice(i) = diag - Covxy.slice(i) + Exy.row(i).t()*Exy.row(i);

    // X - E[X | Y]
    Exy.row(i) = X.row(i) - Exy.row(i);
  }

  // Initial function value

  double F = save_f(B, X, Y, Exy, Covxy, bw, ncore);

  return F;
}


//' @title semi-save solver \code{C++} function
//' @name save_solver
//' @description Sovling the semi-save estimating equations. This is an internal function and should not be called directly.
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param X A matrix of the parameters \code{X}
//' @param Y A matrix of the parameters \code{Y}
//' @param bw Kernel bandwidth for X
//' @param rho (don't change) Parameter for control the linear approximation in line search
//' @param eta (don't change) Factor for decreasing the step size in the backtracking line search
//' @param gamma (don't change) Parameter for updating C by Zhang and Hager (2004)
//' @param tau (don't change) Step size for updating
//' @param epsilon (don't change) Parameter for apprximating numerical gradient, if \code{g} is not given.
//' @param btol (don't change) The \code{$B$} parameter tolerance level
//' @param ftol (don't change) Functional value tolerance level
//' @param gtol (don't change) Gradient tolerance level
//' @param maxitr Maximum number of iterations
//' @param verbose Should information be displayed
//' @references Ma, Y. & Zhu, L. (2012). A semiparametric approach to dimension reduction. Journal of the American Statistical Association, 107(497), 168-179.
//' DOI: \url{https://dx.doi.org/10.1214\%2F12-AOS1072SUPP}.
//' @references Wen, Z. & Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434.
//' DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//'
// [[Rcpp::export]]

List save_solver(arma::mat B,
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

  arma::mat crit(maxitr,3);
  bool invH = true;
  arma::mat eye2P(2*ndr,2*ndr);

  if(ndr < P/2){
    invH = false;
    eye2P.eye();
  }

  // initialize parallel computing

  checkCores(ncore, verbose);

  //precalculate

  arma::mat kernel_matrix_y = KernelDist_multi(Y, ncore, 1);

  arma::rowvec Ky = sum(kernel_matrix_y, 0);

  // X - E[X | Y]
  arma::mat Exy(N, P, arma::fill::zeros);

  // I - cov[X | Y]
  arma::cube Covxy(P, P, N, arma::fill::zeros);
  arma::mat diag = arma::eye(P, P);

#pragma omp parallel for schedule(static) num_threads(ncore)
  for(int i=0; i<N; i++){
    for(int j=0; j<N; j++){
      Exy.row(i) += X.row(j)*kernel_matrix_y(i,j);
      Covxy.slice(i) += (X.row(j).t() * X.row(j)) * kernel_matrix_y(i,j);
    }

    // E[X | Y]
    Exy.row(i) /= Ky(i);

    // E[XX | Y]
    Covxy.slice(i) /= Ky(i);

    // I - cov[X | Y]
    Covxy.slice(i) = diag - Covxy.slice(i) + Exy.row(i).t()*Exy.row(i);

    // X - E[X | Y]
    Exy.row(i) = X.row(i) - Exy.row(i);
  }

  // Initial function value and gradient, prepare for iterations

  double F = save_f(B, X, Y, Exy, Covxy, bw, ncore);

  arma::mat G(P, ndr);
  G.fill(0);
  save_g(B, F, G, X, Y, Exy, Covxy, bw, epsilon, ncore);

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

  if (verbose > 1)
    Rcout << "Initial value,   F = " << F << std::endl;

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

      F = save_f(B, X, Y, Exy, Covxy, bw, ncore);
      save_g(B, F, G, X, Y, Exy, Covxy, bw, epsilon, ncore);

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

    Y_Y = dtX - dtXP;
    SY = std::abs(accu(S % Y_Y));

    if(itr%2 == 0){
      tau = accu(S % S)/SY;
    }else{
      tau = SY/accu(Y_Y % Y_Y);
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
    Rcout << "norm of functional value: " << F << std::endl;
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
