#ifndef ORTHODR_REG_PROJECTION_H
#define ORTHODR_REG_PROJECTION_H

// Called when subtracting already projected coordinates may erase a gap.
static void reg_project_pair(const arma::mat& X, const arma::mat& B,
                             int i, int j, arma::rowvec& delta)
{
  for (arma::uword k = 0; k < B.n_cols; ++k) {
    double sum = 0.0;
    double correction = 0.0;
    for (arma::uword p = 0; p < X.n_cols; ++p) {
      const double a = X(j,p);
      const double b = X(i,p);
      const double diff = a - b;
      const double other = a - diff;
      const double tail = (a - (diff + other)) + (other - b);
      const double product = diff * B(p,k);
      const double next = sum + product;
      const double added = next - sum;
      // Recover roundoff in subtraction, multiplication and summation.
      correction += ((sum - (next - added)) + (product - added))
                  + std::fma(diff, B(p,k), -product) + tail * B(p,k);
      sum = next;
    }
    delta(k) = sum + correction;
  }
}

static arma::mat reg_pair_kernel(const arma::mat& X, const arma::mat& B,
                                 double bw, int ncore)
{
  const int N = X.n_rows;
  const double norm_factor = 1.0 / std::sqrt(static_cast<double>(N) * (N - 1));
  arma::rowvec pair_sd(B.n_cols, arma::fill::zeros);
  arma::rowvec delta(B.n_cols);
  for (int i = 1; i < N; ++i) {
    for (int j = 0; j < i; ++j) {
      reg_project_pair(X, B, i, j, delta);
      // Sum of squared pair distances divided by N(N-1) is the variance.
      for (arma::uword k = 0; k < B.n_cols; ++k)
        pair_sd(k) = std::hypot(pair_sd(k), norm_factor * delta(k));
    }
  }
  const arma::rowvec scale = pair_sd * bw * std::sqrt(2.0);
  arma::mat kernel(N, N, arma::fill::none);
#pragma omp parallel for schedule(guided) num_threads(ncore) if(ncore > 1)
  for (int i = 0; i < N; ++i) {
    kernel(i,i) = 1;
    arma::rowvec distance(B.n_cols);
    for (int j = 0; j < i; ++j) {
      reg_project_pair(X, B, i, j, distance);
      distance /= scale;
      kernel(i,j) = kernel(j,i) = std::exp(-arma::dot(distance, distance));
    }
  }
  return kernel;
}

#endif
