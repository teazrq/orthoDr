# orthoDr (development version)

* Updated direct learning with `intercept = TRUE` to propose adaptive spectral
  step lengths for direction updates. The objective, coefficient penalty,
  mean-dose bounds, initialization, and joint residual convergence criterion
  are unchanged. Fits can still stop before reaching the gradient tolerance;
  inspect `stop.reason` and `gradient.norm`.
* Added optional `intercept = TRUE` for direct dose learning. It learns an
  unpenalized average training dose and penalizes centered kernel coefficients.
  `intercept.bounds` constrains that average, not individual predictions.
  Intercept fits require a positive `lambda.reward` and report their joint
  policy residual. The default remains `intercept = FALSE`; older fitted
  objects retain their prediction behavior.
* Repaired automatic pSAVE initialization for small samples by excluding dose
  groups without a usable predictor covariance and handling the dependency's
  one-breakpoint slicing case. Structurally uninformative slice configurations
  now give a specific error suggesting `B.initial`. Valid historical thresholds,
  including ties, retain their contributions for both dose-learning methods.
* Revised direct dose learning to optimize one penalized reward objective over
  both dose-rule coefficients and directions. `lambda` now controls its ridge
  initialization only; the new `lambda.reward` controls the squared coefficient
  penalty. Kernel centres and projected scales follow the fitted directions.
  Direct fits now report the estimated reward, penalty, both gradient norms,
  and stopping reason. Failed line searches preserve the last completed sweep,
  and small steps no longer imply convergence. The direct `fn` field includes
  the reward penalty, so it cannot be compared directly with older values.
  Existing fitted objects retain their prediction behavior.
* Expanded the tutorials with reduced predictors, predicted survival curves,
  and explanations of the different targets and dimensions in dose finding.
  Clarified prediction inputs and corrected example and reference descriptions.
* Corrected the time-axis labels in 3D survival plots to match the transformed
  time coordinates. The plotted survival probabilities are unchanged.
* Added a pkgdown GitHub Pages site with installation instructions, four
  workflow tutorials, the function reference, changelog, and checked example
  downloads. Site builds reuse the tested guide and deploy through GitHub Actions.
* Added a versioned getting-started guide with runnable regression, survival,
  dose, prediction, and optimizer examples, exercised during package checks.
  Corrected the four legacy plotting calls to use `FUN` and retained their
  script in the repository.
* Corrected regression and survival prediction coordinates. New forward
  survival fits retain bandwidths; older forward fits require the original
  bandwidth or refitting before prediction.
* Corrected survival product-limit risk denominators and time labels, and
  made tied-time risk sets consistent across initialization and all three
  estimating equations. Gaussian weight underflow now produces a warning
  and an undefined curve instead of an all-one curve.
* Used a common observed-dose scale for dose kernels, resolved the default
  method, and made fractional and single-point dose grids consistent.
* Protected optimizer step calculations at stationary and degenerate updates,
  prevented numerical stagnation from reporting movement-based convergence,
  and corrected completed iteration counts and maximization histories.
  Direct-dose convergence now requires stationarity in both parameter blocks.
* Validated general optimizer callback values on every evaluation and accepted
  integer as well as double numeric callbacks.
* Preserved requested regression initializer dimensions, repaired CP-SIR and
  survival matrix shapes for one feature or event, and fixed pSAVE indexing,
  documented hMAVE vector inputs, and distance argument handling.
* These numerical corrections can change fitted results and predictions.
  Analyses affected by these issues should be rerun.

* Added `testthat` test suite for the general optimizer (`ortho_optim`)
* Fixed CRLF line endings in source files for cross-platform compatibility
* Fixed const-reference issues in C++ headers
* Migrated to roxygen2-managed NAMESPACE with explicit exports
* Updated CITATION to modern `bibentry()` API

# orthoDr 0.6.2

* Added personalized dose finding models.

# orthoDr 0.5.1

* Fixed minor compiling issues.

# orthoDr 0.5.0

* Added additional regression functions.
* Added additional features for kernel prediction based on the dimensionality reduced space.

# orthoDr 0.3.0

* Added parallel computing using OpenMP.

# orthoDr 0.2.0

* Improved speed.

# orthoDr 0.1.0

* First release on CRAN.
