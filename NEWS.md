# orthoDr (development version)

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
  Direct-dose gradient stopping retains its conditional-policy interpretation.
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

# orthoDr 1.0.0

* CRAN release v1.0

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
