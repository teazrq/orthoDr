#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

double gen_f(arma::mat &B, Rcpp::Function f, Environment env)
{
  SEXP x = Rcpp_eval(f(B), env);
  return(REAL(x)[0]);
}


void gen_g(arma::mat B, arma::mat &G, Rcpp::Function g, Environment env)
{
  SEXP x = Rcpp_eval(g(B), env);

  int P = B.n_rows;
  int ndr = B.n_cols;

  for (int k=0; k<ndr; k++)
    for (int j=0; j<P; j++)
      G(j, k) = REAL(x)[j + k*P];

  return;
}


void gen_g_approx(arma::mat &B, arma::mat &G, Rcpp::Function f, Rcpp::Function g, Environment env, double epsilon)
{

  double F0 = gen_f(B, f, env);
  int P = B.n_rows;
  int ndr = B.n_cols;
  double temp;

  for (int j = 0; j < ndr; j++)
  {
    for(int i = 0; i < P; i++)
    {
      // small increment
      temp = B(i,j);
      B(i,j) += epsilon;

      // calculate gradient
      G(i,j) = (gen_f(B, f, env) - F0) / epsilon;

      // reset
      B(i,j) = temp;

    }
  }

  return;
}

double dmax(double a, double b);
double dmin(double a, double b);


//' @title General solver \code{C++} function
//' @name gen_solver
//' @description A general purpose optimization solver with orthogonality constraint. For details, please see the original \code{MATLAB} code by Wen and Yin (2013). This is an internal function and should not be called directly.
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param f A function that calculates the objective function value. The first argument should be \code{B}. Returns a single value.
//' @param g A function that calculates the gradient. The first argument should be \code{B}. Returns a matrix with the same dimension as \code{B}. If not specified, then the numerical approximation is used.
//' @param env Environment passed to the Rcpp function for evaluating \code{f} and \code{g}
//' @param useg If true, the gradient is calculated using \code{g} function, otherwise numerically approximated
//' @param rho (don't change) Parameter for control the linear approximation in line search
//' @param eta (don't change) Factor for decreasing the step size in the backtracking line search
//' @param gamma (don't change) Parameter for updating C by Zhang and Hager (2004)
//' @param tau (don't change) Step size for updating
//' @param epsilon (don't change) Parameter for approximating numerical gradient, if \code{g} is not given.
//' @param btol (don't change) The \code{$B$} parameter tolerance level
//' @param ftol (don't change) Functional value tolerance level
//' @param gtol (don't change) Gradient tolerance level
//' @param maxitr Maximum number of iterations
//' @param verbose Should information be displayed
//' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//' @references Zhang, H. and Hager, W. W., "A nonmonotone line search technique and its application to unconstrained optimization." SIAM J. Optim. 14 (2004): 1043–1056. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//' @examples
//' # This function should be called internally. When having all objects pre-computed, one can call
//' # gen_solver(B, f, g, env, useg, rho, eta, gamma, tau, epsilon, btol, ftol, gtol, maxitr, verbose)
//' # to solve for the parameters B.
// [[Rcpp::export]]


List gen_solver(arma::mat B,
                 Rcpp::Function f,
                 Rcpp::Function g,
                 Environment env,
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
                 int verbose)
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

  // Initial function value and gradient, prepare for iterations

  double F = gen_f(B, f, env);

  arma::mat G(P, ndr);
  G.fill(0);

  if (useg)
    gen_g(B, G, g, env);
  else
    gen_g_approx(B, G, f, g, env, epsilon);

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
  arma::mat Y;
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

    // line search
    while(true){
      if(invH){
        diag_n.eye();
        B = solve(diag_n + tau * H, BP - tau * RX);
      }else{
        aa = solve(eye2P + 0.5 * tau * VU, VX);
        B = BP - U * (tau * aa);
      }

      F = gen_f(B, f, env);

      if (useg)
        gen_g(B, G, g, env);
      else
        gen_g_approx(B, G, f, g, env, epsilon);

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
    }else{ // algorithm 1
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

    Y = dtX - dtXP;
    SY = std::abs(accu(S % Y));

    if(itr%2 == 0){
      tau = accu(S % S)/SY;
    }else{
      tau = SY/accu(Y % Y);
    }

    tau = dmax(dmin(tau, 1e20), 1e-20);
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
