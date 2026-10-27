# orthoDr

[![CRAN Version](https://www.r-pkg.org/badges/version/orthoDr)](https://CRAN.R-project.org/package=orthoDr)

orthoDr estimates a small set of linear combinations of predictors for
regression, survival analysis, and personalized dose finding. Use the fitted
directions to obtain reduced predictors, predict responses or survival
probabilities, or estimate a dose rule. The package also provides a general
optimizer for problems with orthogonality constraints.

## Features

- **Regression**: SIR, SAVE, PHD, locally efficient, and efficient estimation
- **Survival**: dimension reduction for censored outcomes
- **Dose finding**: personalized dose estimation via partial SAVE and direct/pseudo-direct learning
- **General optimizer**: `ortho_optim()` for any orthogonality-constrained objective

## Installation

This website documents **development version 2.0.0**. Install it from GitHub
to follow the tutorials:

```r
install.packages("remotes") # once, if needed
remotes::install_github("teazrq/orthoDr")
```

Installation from GitHub requires R's C++ build tools. See
[Get Started](https://teazrq.github.io/orthoDr/articles/get-started.html) for details.
For the CRAN release, version 0.6.8, use `install.packages("orthoDr")` and its
[reference manual](https://cran.r-project.org/web/packages/orthoDr/refman/orthoDr.html).
Numerical results can differ between versions.

## Getting started

- [Get Started](https://teazrq.github.io/orthoDr/articles/get-started.html): installation
  and tutorials for regression, survival, personalized dose finding, and the general optimizer.
- [Runnable guide](https://teazrq.github.io/orthoDr/examples/getting-started.R) and
  [supplemental survival example](https://teazrq.github.io/orthoDr/examples/survival-optimizer.R),
  including the bundled melanoma data and plots of reduced predictors.
- [Development function reference](https://teazrq.github.io/orthoDr/reference/index.html),
  [release notes](https://teazrq.github.io/orthoDr/news/index.html), and
  [CRAN 0.6.8 reference manual](https://cran.r-project.org/web/packages/orthoDr/refman/orthoDr.html).

For documentation matching your installed version, use `help(package = "orthoDr")`.

For direct dose learning, `lambda` controls ridge initialization and
`lambda.reward` controls the squared coefficient penalty when optimizing
estimated reward. Inspect `converge`, `gradient.norm`, and `stop.reason`
before using a fitted rule. The nominal reward penalty depends on dose units
and reward scale; it is not selected by the initialization's GCV procedure.
Use `intercept = TRUE` to learn an unpenalized average training dose while
penalizing the centered kernel coefficients. Its `intercept.bounds` constrain
that average, not individual predictions. This option requires a positive
`lambda.reward`; the default remains `intercept = FALSE`.

## References

- Wen, Z. & Yin, W. (2013). A feasible method for optimization with orthogonality constraints. *Mathematical Programming*, 142(1-2), 397–434.
- Ma, Y. & Zhu, L. (2012). A semiparametric approach to dimension reduction. *JASA*, 107(497), 168–179.
- Ma, Y. & Zhu, L. (2013). Efficient estimation in sufficient dimension reduction. *Annals of Statistics*, 41(1), 250–268.
- Sun, Q., Zhu, R., Wang, T. & Zeng, D. (2017). Counting process based dimension reduction for censored outcomes. arXiv:1704.05046.
- Zhou, W., Zhu, R. & Zeng, D. (2021). A parsimonious personalized dose-finding model via dimension reduction. *Biometrika*, 108(3), 643-659. [doi:10.1093/biomet/asaa087](https://doi.org/10.1093/biomet/asaa087).
