#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]


double surv_f(arma::mat B,
              arma::mat X,
              arma::mat Phit,
              NumericMatrix inRisk,
              double bw,
              NumericVector Fail_Ind)
{
  // This function computes the estimation equations and its 2-norm for the survival dimensional reduction model
  // It only implement the dN method, with phi(t)

  int N = X.n_rows;
  int P = X.n_cols;
  int nFail = inRisk.ncol();

  NumericMatrix BX_NM = wrap(X * B);

  NumericMatrix kernel_matrix(N, N);

  for (int i = 0; i < N; i++)
  {
    kernel_matrix(i, i) = 1;
    for (int j = i + 1; j < N; j ++)
    {
      kernel_matrix(j,i) = exp(-sum(pow(BX_NM.row(i)-BX_NM.row(j),2))/bw/bw);
      kernel_matrix(i,j) = kernel_matrix(j,i);
    }
  }

  NumericMatrix TheCenter(nFail, P);
  NumericMatrix X_NM = wrap(X);
  NumericVector TheCond(P);

  double unweighted_sum;
  double weights;

  for(int j=0; j<nFail; j++)
  {
    for(int i=0; i<P; i++)
    {
      unweighted_sum = 0;
      weights = 0;

      for(int k=0; k<N; k++)
      {
        if(inRisk(k,j)==true)
        {
          unweighted_sum = unweighted_sum + X_NM(k,i) * kernel_matrix(k,Fail_Ind[j]-1);
          weights = weights + kernel_matrix(k,Fail_Ind[j]-1);
        }
      }
      TheCond[i] = unweighted_sum/weights;
    }
    TheCenter(j,_) = X_NM(Fail_Ind[j]-1,_) - TheCond;
  }

  arma::mat matrixsum(P, P);
  matrixsum.fill(0);
  arma::mat TheCenter_arma = as<arma::mat>(TheCenter);

  for(int i=0; i<nFail; i++)
  {
    matrixsum = matrixsum + Phit.col(i) * TheCenter_arma.row(i);
  }

  NumericMatrix matrixsum_NM = wrap(matrixsum);
  double ret = sum(pow(matrixsum_NM/N,2));

  return ret;
}


void surv_g(arma::mat B,
            arma::mat& G,
            arma::mat X,
            arma::mat Phit,
            NumericMatrix inRisk,
            double bw,
            NumericVector Fail_Ind,
            double epsilon)
{
  // This function computes the gradiant of the estimation equations

  arma::mat B_new = B;
  double F0 = surv_f(B, X, Phit, inRisk, bw, Fail_Ind);
  int P = B.n_rows;
  int ndr = B.n_cols;

  for (int j = 0; j < ndr; j++)
  {
    for(int i = 0; i < P; i++)
    {
      // small increment
      B_new(i,j) = B(i,j) + epsilon;

      // calculate gradiant
      G(i,j) = (surv_f(B_new, X, Phit, inRisk, bw, Fail_Ind) - F0) / epsilon;

      // reset
      B_new(i,j) = B(i,j);
    }
  }

  return;
}


double dmax(double a, double b)
{
  if (a > b)
    return a;

  return b;
}

double dmin(double a, double b)
{
  if (a > b)
    return b;

  return a;
}


//' @title surv_solver
//' @name surv_solver
//' @description The main optimization function for survival dimensional reduction, the IR-CP method. This is an internal function and should not be called directly.
//' @keywords internal
//' @param B A matrix of the parameters \code{B}, the columns are subject to the orthogonality constraint
//' @param X The covariate matrix
//' @param Phit Phit as defined in Sun et al. (2017)
//' @param inRisk A matrix of indicators shows whether each subject is still alive at each time point
//' @param bw A Kernel bandwidth, assuming each variable have unit variance
//' @param Fail_Ind The locations of the failure subjects
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
//' @references Sun, Q., Zhu, R., Wang T. and Zeng D. "Counting Process Based Dimension Reduction Method for Censored Outcomes." (2017) \url{https://arxiv.org/abs/1704.05046} .
//' @references Wen, Z. and Yin, W., "A feasible method for optimization with orthogonality constraints." Mathematical Programming 142.1-2 (2013): 397-434. DOI: \url{https://doi.org/10.1007/s10107-012-0584-1}
//' @examples
//' # This function should be called internally. When having all objects pre-computed, one can call
//' # surv_solver(B, X, Phit, inRisk, bw, Fail.Ind,
//' #             rho, eta, gamma, tau, epsilon, btol, ftol, gtol, maxitr, verbose)
//' # to solve for the parameters B.
//'
// [[Rcpp::export]]

List surv_solver(arma::mat B,
        							arma::mat X,
        							arma::mat Phit,
        							NumericMatrix inRisk,
        							double bw,
        							NumericVector Fail_Ind,
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

  int ndr = B.n_rows;
  int P = B.n_cols;

  NumericMatrix crit(maxitr,3);
  bool invH = true;
  arma::mat eye2P(2*P,2*P);

  if(P < ndr/2){
    invH = false;
    eye2P.eye();
  }

  // Initial function value and gradient, prepare for iterations

  double F = surv_f(B, X, Phit, inRisk, bw, Fail_Ind);

  arma::mat G(ndr, P);
  G.fill(0);
  surv_g(B, G, X, Phit, inRisk, bw, Fail_Ind, epsilon);

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
  arma::mat diag_n(ndr, ndr);
  arma::mat aa;
  arma::mat S;
  double BDiff;
  double FDiff;
  arma::mat Y;
  NumericMatrix S_Y;
  double SY;
  NumericMatrix S_S;
  NumericMatrix Y_Y;

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

      F = surv_f(B, X, Phit, inRisk, bw, Fail_Ind);
      surv_g(B, G, X, Phit, inRisk, bw, Fail_Ind, epsilon);

      if (verbose > 1 && (itr % 10 == 0) )
        Rcout << "At iteration " << itr << ", F = " << F << std::endl;

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
    BDiff = norm(S, "fro")/sqrt((double) ndr);
    FDiff = std::abs(FP - F)/(std::abs(FP)+1);

    Y = dtX - dtXP;
    S_Y = wrap(S % Y);
    SY = std::abs(sum(S_Y));

    if(itr%2 == 0){
      S_S = wrap(S % S);
      tau = sum(S_S)/SY;
    }else{
      Y_Y = wrap(Y % Y);
      tau = SY/sum(Y_Y);
    }

    tau = dmax(dmin(tau, 1e10), 1e-20);
    crit(itr-1,0) = nrmG;
    crit(itr-1,1) = BDiff;
    crit(itr-1,2) = FDiff;


    if (itr >= 5) // so I will run at least 5 iterations before checking for convergence
    {
      NumericMatrix mcrit(5, 3);
      for (int i=0; i<5; i++)
      {
        mcrit(i, _) = crit(itr-i,_);
      }

      if ( (BDiff < btol && FDiff < ftol) || (nrmG < gtol) || ((mean(mcrit(_,1)) < btol) && (mean(mcrit(_,2)) < ftol)) )
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

  arma::mat diag_P(P,P);
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
